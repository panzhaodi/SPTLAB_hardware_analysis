import pylab as pl
import numpy as np
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from glob import glob
import sptpol_software.util.time as sptTime
#'/home/cryo/Data/output/20140414_cooldown/20140329_214118_IV'
def plotIv( given ,channel,resistance, parRes,figure=1,note='',
           plot_iv = False,
           plot_loopgain = False,
           plot_by_squid = True,
           plot_by_bolo = True,
           return_data = False):
 
    item = dict()
    # figure out what type of object 'given' is and process it
     
    pkls = glob('/home/cryo/Data/output/20140511_cooldown/temp/'+given+'/*TuneBoloCombDan*.pkl')
    if pkls:
            for _pkl in pkls:
                _label = _pkl.rpartition('/')[-1].partition('_')[0]
                item[_label] = [pkl.load(open(_pkl,'r'))]
     
    
    # Setup some specialized subplots 
    num_subplots = 1
    # Some features are not compatible with plot_by_bolo, 
    #    turn off these features here
  
    # no matter what we will be doing we will use figure(figure), 
    # so we will initiate it and clear it now
    pl.figure(figure)
    pl.clf()
    
    data = dict()
    max_n_sq = len(item)
    # this loop is over all pkls (1 pkl = 1 squid, so looping over all squids)
    
    for _n_sq,_sq in enumerate(sorted(item.keys())):
    
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
                now_subplot = 1
                v_subplot = 1
                h_subplot = 1
                pl.subplot(v_subplot,h_subplot,now_subplot)
                RP_plotnum = now_subplot
                # label the axes
                pl.xlabel('Power [pW]')
                pl.ylabel('Resistance [Ohm]')
                pl.grid()
                pl.title(_sq.partition('-')[0]+'_'+'Channel'+channel)

            else:
                v_subplot = 1
                h_subplot = now_subplot
                pl.subplot(v_subplot,h_subplot,now_subplot)
                RP_plotnum = now_subplot
            # actually plot hte RP curve
            if str(chan)==channel:
                pl.plot(vals['p'],vals['r'],'.-',label=_sq.rpartition('-')[-1])
                pl.legend(loc='best')
                print 'The power of '+_sq.partition('-')[0]+' Channel '+channel+ ' at '+_sq.rpartition('-')[-1]+' %f Ohm is %f pW'  %(resistance, np.interp(resistance, vals['r'], vals['p']))   
               # print( 'The power of'+_sq.partition('-')[0]++'Channel'+channel+'at'+ label+ 'and'+'%f Ohm is %f MuW'  %(resistance, np.interp(resistance, vals['r'], vals['p']))   )
        # label plots if plot_by_bolo = False
         
        pl.tight_layout()

   



