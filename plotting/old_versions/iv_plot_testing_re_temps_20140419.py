import pylab as pl
import numpy as np
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from glob import glob
import re
import sptpol_software.util.time as sptTime

#cold_template = 'Tcoldload = (\d+)'
cold_template = 'Tcoldload = *'

def plotIv(given, figure=1,note='',parRes=0.03,
           plot_iv = False,
           plot_loopgain = False,
           plot_by_squid = False,
           plot_by_bolo = True,
           return_data = False):
    '''
    Quickly plot bolometer IV data

    INPUTS:
       given :(dict, directory, .pkl location) 
              dict - all iv data from the dict will be plotted (formatting must 
                  be correct). 
              directory - All IV pkl files in the directory will be read in and 
                  plotted. 
              .pkl location - The pkl will be readin and the IV data it contains 
                  will be plotted. 
                ************************** Currently feeding a pkl location is broken ******************************

       figure [1]: (int) The figure to plot the data on, if the figure already 
              exists it will be cleared. If plot_by_squid=True then figures from
              [figure : figure+n_squids] will be used. 

       parRes [0.03] (float) The parasitic resistance to use for plotting up the 
              IV data. 

       plot_iv [True]: (bool) Include an iv subplot in raw readout units. 

       plot_loopgain [True]: (bool) Include a subplot of the loopgain derivied
              from the IV data.
        
       plot_by_squid [False]: (bool) Creates a plot for each of the squids. 
       
       plot_by_bolo [False]: (bool) Creates a subplot of PvsR for each bolo within 
              a squid plot. This option overrides the following commands;
                 plot_by_squid = True
                 plot_iv = False
                 plot_loopgain = False

       return_data [False]: (bool) Return a dict including the IV data in actual units.

    OUTPUT:
       A dict containing IV data coverted to actual units (only if return_data=True).
    '''

    item = dict()
    # figure out what type of object 'given' is and process it
    try:
        given.keys()
        item[''] = given
    except AttributeError:
        pkls = glob(given+'/*TuneBoloCombDan*.pkl')
        if pkls:
            for _pkl in pkls:
                _label = _pkl.rpartition('/')[-1].partition('_')[0]
                item[_label] = [pkl.load(open(_pkl,'r'))]
            # if a cycle_tune.log file is present, extract the temp of the coldload 
            #    and bolos from it
            logfile = glob(given+'/cycle_tune.log')[0]
            if logfile:
                logfile = open(logfile,'r')
                Tbolos_all,Tcoldload_all = [],[]
                for ln in logfile:
                    Tbolos_now = re.search(r'bolos.=.*',ln)
                    Tcoldload_now = re.search(r'ColdLoad.=.*',ln)
                    if Tcoldload_now:
                        Tcoldload_all.append(float(Tcoldload_now.group().rpartition(' ')[-1]))
                    if Tbolos_now:
                        Tbolos_all.append(float(Tbolos_now.group().rpartition(' ')[-1]))
                Tcoldload,Tcoldload_err = np.average(Tcoldload_all), (max(Tcoldload_all)-min(Tcoldload_all))/2
                Tbolos,Tbolos_err = np.average(Tbolos_all), (max(Tbolos_all)-min(Tbolos_all))/2
                print 'Tcoldload = ',Tcoldload,' +/- ',Tcoldload_err 
                print 'Tbolos = ',Tbolos,' +/- ',Tbolos_err 
                
        else:
            print 'reading pickle'
            # if it does not have keys or pkls in it, 
            #   assume it is a pkl file and read it in
            item = [pkl.load(open(given,'r'))]
    
    # Setup some specialized subplots 
    num_subplots = 1
    if plot_iv: num_subplots+=1
    if plot_loopgain: num_subplots+=1

    # Some features are not compatible with plot_by_bolo, 
    #    turn off these features here
    if plot_by_bolo:
        plot_by_squid = True
        plot_iv = False
        plot_loopgain = False

    # no matter what we will be doing we will use figure(figure), 
    # so we will initiate it and clear it now
    pl.figure(figure)
    pl.clf()
    
    data = dict()
    max_n_sq = len(item)
    # this loop is over all pkls (1 pkl = 1 squid, so looping over all squids)
    for _n_sq,_sq in enumerate(sorted(item.keys())):
        if plot_by_squid:
            pl.figure(figure+_n_sq)
            pl.clf()
        _it = item[_sq][0]
        data[_sq] = dict()
        # this loop is over all bolos within a squid
        for _n_chan, chan in enumerate(sorted(_it['data'].keys())):
            # convert readout units to physical units
            vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'
            data[_sq]['Ch'+str(chan)] = dict({'data':_it['data'][chan],
                                              'vals':vals})
            now_subplot = 1
            # there is a lot of messy subplot bookkeeping below
            if plot_iv:
                # add the iv plot and fill it
                pl.subplot(1,num_subplots,now_subplot)
                pl.plot(vals['v'],vals['i'],label=str(chan)+''+'')
                # pl.plot(vals['vRaw'],vals['iRaw'],label=str(chan)+'_raw')
                # uncomment the line above if you are interested in plotting uncorrected 
                #    IV curves also
                iv_plotnum = now_subplot
                now_subplot+=1
            if plot_loopgain:
                # add the loopgain subplot and fill it
                pl.subplot(1,num_subplots,now_subplot)
                pl.plot(vals['p'],vals['loopgain'],label=str(chan)+''+'')
                loopgain_plotnum = now_subplot
                now_subplot+=1

            if plot_by_bolo: # plot_by_bolo will skip right to here for plotting
                now_subplot = _n_chan+1
                v_subplot = 2
                h_subplot = 4
                pl.subplot(v_subplot,h_subplot,now_subplot)
                RP_plotnum = now_subplot
                # label the axes
                pl.xlabel('Power [pW]')
                pl.ylabel('Resistance [Ohm]')
                pl.grid()
                pl.title(_sq+'_'+str(chan))

            else:
                v_subplot = 1
                h_subplot = now_subplot
                pl.subplot(v_subplot,h_subplot,now_subplot)
                RP_plotnum = now_subplot
            # actually plot the RP curve
            if logfile:
                RP_label = _sq+'_'+str(chan)+' %.3gK' %Tcoldload
            else:
                RP_label = _sq+'_'+str(chan)
            pl.plot(vals['p'],vals['r'],'.-',label=RP_label)

        # label plots if plot_by_bolo = False
        if plot_by_squid or _n_sq+1 == max_n_sq:
            if plot_iv:
                pl.subplot(v_subplot,h_subplot,iv_plotnum)
                pl.xlabel('Voltage [uV]')
                pl.ylabel('Current [uI]')
                pl.legend(loc=2)
                pl.grid()
            if plot_loopgain:
                pl.subplot(v_subplot,h_subplot,loopgain_plotnum)
                pl.xlabel('Power [pW]')
                pl.ylabel('loopgain [uV]')
                pl.legend(loc=2)
                pl.grid()
            if not plot_by_bolo:
                pl.subplot(v_subplot,h_subplot,RP_plotnum)
                pl.xlabel('Power [pW]')
                pl.ylabel('Resistance [Ohm]')
                pl.grid()
                pl.legend(loc = 4)

        pl.tight_layout()

    if return_data:
        return data


def plotIvEtf(given, figure=1,note='',parRes=0.03,
              ownPlots=False,plot_iv=True,
              returnData=False,smooth_loop=False):

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
            if smooth_loop:
                num_loop_points_to_smooth = 5
                if _frac < 0.99:
                    num_loop_points_to_smooth = 7
                poly_vals= np.polyfit(vals['p'][:num_loop_points_to_smooth],
                                      vals['r'][:num_loop_points_to_smooth],2)
#                vals['smooth_loopgain'] = np.polyfit(vals['p'],vals['loopgain'],2)
                poly_now = np.poly1d(poly_vals)
                tight_power = np.linspace(min(vals['p'][:num_loop_points_to_smooth]),max(vals['p'][:num_loop_points_to_smooth]),200)
                tight_resist = poly_now(tight_power)
                pl.plot(tight_power,tight_resist,color = pl.gca().lines[-1].get_color(),alpha=0.7)
                tight_loopgain = (tight_power[:-1]/tight_resist[:-1])*(np.diff(tight_resist)/np.diff(tight_power))

            ax = pl.subplot(133)
            pl.plot(vals['p'],vals['loopgain'],'.-',alpha=0.3,#label=str(_frac),
                    color = pl.cm.gist_rainbow(float(rfrac_count-1-_nfrac)/rfrac_count))
            if smooth_loop:
                pl.plot(vals['p'][0],vals['loopgain'][0],'o',label=str(_frac)+' : Loop=%.3g : %.3g'
                        %(vals['loopgain'][0],tight_loopgain[0]),
                    color = pl.gca().lines[-1].get_color())
                vals['final_smooth_loopgain'] = tight_loopgain[0]
            else:
                pl.plot(vals['p'][0],vals['loopgain'][0],'o',label=str(_frac)+' : Loop=%.3g'%vals['loopgain'][0],
                    color = pl.gca().lines[-1].get_color())

            if smooth_loop:
                pl.plot(tight_power[:-1],tight_loopgain,
                        color = pl.gca().lines[-1].get_color())
                pl.plot(tight_power[0],tight_loopgain[0],'*',
                        color = pl.gca().lines[-1].get_color())
                loop2 = tight_loopgain[0]

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
