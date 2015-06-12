import pylab as pl
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from glob import glob
import sptpol_software.util.time as sptTime

def plotIv(given, figure=1,note='',parRes=0.03,
           plot_iv=True,
           plot_loopgain=True,
           return_data=False):

    item = dict()
    try:
        given.keys()
        item[''] = given
    except AttributeError:
        pkls = glob(given+'/*TuneBoloCombDan*.pkl')
        if pkls:
            for _pkl in pkls:
                _label = _pkl.rpartition('/')[-1].partition('_')[0]
                item[_label] = [pkl.load(open(_pkl,'r'))]
        else:
            print 'reading pickle'
            # if it does not have keys or pkls in it, 
            #   assume it is a pkl file and read it in
            item = [pkl.load(open(item,'r'))]
    
    num_subplots = 1
    if plot_iv: num_subplots+=1
    if plot_loopgain: num_subplots+=1
 
    pl.figure(figure)
    pl.clf()
    data = dict()
    for _sq in item:
        _it = item[_sq][0]
        data[_sq] = dict()
        for chan in _it['data'].keys():
            vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'
            data[_sq]['Ch'+str(chan)] = dict({'data':_it['data'][chan],
                                              'vals':vals})
            now_subplot = 1
            if plot_iv:
                pl.subplot(1,num_subplots,now_subplot)
                pl.plot(vals['v'],vals['i'],label=str(chan)+''+'')
                pl.plot(vals['vRaw'],vals['iRaw'],label=str(chan)+'_raw')
                iv_plotnum = now_subplot
                now_subplot+=1
            if plot_loopgain:
                pl.subplot(1,num_subplots,now_subplot)
                pl.plot(vals['p'],vals['loopgain'],label=str(chan)+''+'')
                loopgain_plotnum = now_subplot
                now_subplot+=1
            pl.subplot(1,num_subplots,now_subplot)
            pl.plot(vals['p'],vals['r'],'.-',label=_sq+'_'+str(chan))
            RP_plotnum = now_subplot
#        pl.plot(vals['vRaw'],vals['iRaw'],label=str(chan)+'_raw')

#        pl.suptitle('BoloIV -- Channel'+str(chan)+' '+note)
    if plot_iv:
        pl.subplot(1,num_subplots,iv_plotnum)
        pl.xlabel('Voltage [uV]')
        pl.ylabel('Current [uI]')
        pl.legend(loc=2)
        pl.grid()
        
    if plot_loopgain:
        pl.subplot(1,num_subplots,loopgain_plotnum)
        pl.xlabel('Power [pW]')
        pl.ylabel('loopgain [uV]')
        pl.legend(loc=2)
        pl.grid()

    pl.subplot(1,num_subplots,RP_plotnum)
    pl.xlabel('Power [pW]')
    pl.ylabel('Resistance [Ohm]')
    pl.grid()
    pl.legend(loc = 4)
    pl.tight_layout()

    if return_data:
        return data


def plotIvEtf(given, figure=1,note='',parRes=0.03,
              ownPlots=False,plot_iv=True,
              returnData=False):

    ''' feed me a directory that took etf data
    '''
    
    pkls = glob(given+'/*TuneBoloCombDan*.pkl')
    if not pkls:
        print 'There are no TuneBoloCombDan pkls in:'+given
        return
    # get a list of all the squids
    squids = list(set([_pkl.rpartition('/')[-1].rpartition('_Tune')[0] for _pkl in pkls]))

    item = dict()
    for _sq in squids:
        item[_sq] = dict()
    for _pkl in pkls:
        _sq = _pkl.rpartition('/')[-1].partition('_')[0]
        _label = sptTime.toSptDatetime(_pkl.rpartition('_')[-1].rpartition('.pkl')[0]).fil
        item[_sq][_label] = [pkl.load(open(_pkl,'r'))]
        
    bolos = []
    data = dict()
    outDict = dict()

    pl.figure(figure)
    pl.clf()
    for _sq in item:
        for _file in item[_sq]:
            _now = item[_sq][_file][0]
            _frac = _now['arguments']['channel_frac_resistance']
            for _i, _chan in enumerate(_now['data'].keys()):
                try:
                    complete_crap = data[_sq+'Ch'+str(_chan)]
                except KeyError:
                    data[_sq+'Ch'+str(_chan)] = dict()
                
                vals = iv.analyzeIVData(_now['data'][_chan],iv_type='TuneBoloCombDan', 
                                        parRes = parRes, cgain = _now['carrier_gain'], ngain = _now['nuller_gain'])
                data[_sq+'Ch'+str(_chan)][_frac[_i]] = vals


    for _i,_bolo in enumerate(data.keys()):
        pl.figure(_i+figure);pl.clf()
        pl.suptitle(_bolo+' Rpar='+str(parRes))
        outDict[_bolo] = dict({'Rparasitic':parRes,'Rfrac':dict()})
        rfrac_count = len(data[_bolo].keys())
        for _nfrac,_frac in enumerate(sorted(data[_bolo].keys())):
            vals = data[_bolo][_frac]
            outDict[_bolo]['Rfrac'][_frac] = vals
            ax = pl.subplot(131)
            pl.plot(vals['v'],vals['i'],alpha=0.3,#label=str(_frac),
                    color = pl.cm.gist_rainbow(float(rfrac_count-1-_nfrac)/rfrac_count))
            pl.plot(vals['v'][0],vals['i'][0],'o',label=str(_frac),
                    color = pl.gca().lines[-1].get_color())
#                pl.plot(vals['vRaw'],vals['iRaw'],label=str(frac))
            ax = pl.subplot(132)
            pl.plot(vals['p'],vals['r'],'.-',alpha=0.3,#label=str(_frac),
                    color = pl.cm.gist_rainbow(float(rfrac_count-1-_nfrac)/rfrac_count))
            pl.plot(vals['p'][0],vals['r'][0],'o',label=str(_frac),
                    color = pl.gca().lines[-1].get_color())
            ax = pl.subplot(133)
            pl.plot(vals['p'],vals['loopgain'],'.-',alpha=0.3,#label=str(_frac),
                    color = pl.cm.gist_rainbow(float(rfrac_count-1-_nfrac)/rfrac_count))
            pl.plot(vals['p'][0],vals['loopgain'][0],'o',label=str(_frac)+' : Loop=%.3g'%vals['loopgain'][0],
                    color = pl.gca().lines[-1].get_color())

        ax = pl.subplot(131)
        pl.xlabel('Voltage [uV]')
        pl.ylabel('Current [uI]')
        hands,llabels = ax.get_legend_handles_labels()
        pl.legend(hands[::-1],llabels[::-1],loc=2)
        pl.grid()
        ax = pl.subplot(132)
        pl.xlabel('Power [pW]')
        pl.ylabel('Resistance [Ohm]')
        hands,llabels = ax.get_legend_handles_labels()
        pl.legend(hands[::-1],llabels[::-1],loc=4)
        pl.grid()
        ax = pl.subplot(133)
        pl.xlabel('Power [pW]')
        pl.ylabel('Loopgain [P/R * dR/dP]')
        hands,llabels = ax.get_legend_handles_labels()
        pl.legend(hands[::-1],llabels[::-1],loc=1)
        pl.grid()
        pl.tight_layout()
    if returnData:
        return outDict
