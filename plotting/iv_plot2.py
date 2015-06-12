import pylab as pl
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from glob import glob

def get_iv(sq,chan,temp='???',parRes=0.03,plot_iv=True,legends=False):
	"Example call: get_iv('Sq4',2,'300')"

	path = '/home/cryo/Data/output/20131220_cooldown/br103/temp_plot'
	item=dict()

	pkls = glob(path+'/*'+temp+'mK_IV/'+sq+'_TuneBoloCombDan*.pkl')
	if len(pkls)==0:
	    print "No matching .pkl file found."
	    return
	pkls.sort()	
	for _pkl in pkls:
	    _label= _pkl.partition('mK')[0].rpartition('/')[-1]
            item[_label] = [pkl.load(open(_pkl,'r'))]	

	pl.figure()
    	pl.clf()

	for _run in item:
            _it = item[_run][0]
            vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

	    if plot_iv:
                pl.subplot(121)
                pl.plot(vals['v'],vals['i'],label=_run.rpartition('_')[-1]+'mK')
                #pl.plot(vals['vRaw'],vals['iRaw'],label=_run.rpartition('_')[-1]+'mK_raw')
                pl.subplot(122)
                pl.plot(vals['p'],vals['r'],'.-',label=_run.rpartition('_')[-1]+'mK')
	
	
	pl.subplot(121)
	pl.xlabel('Voltage [uV]')
        pl.ylabel('Current [uI]')
	if legends:        
	    pl.legend()
	pl.subplot(122)
    	pl.xlabel('Power [pW]')
    	pl.ylabel('Resistance [Ohm]')
	if legends:    	
	    pl.legend(loc = 4)
	pl.suptitle(sq +'_Ch'+str(chan))

	pl.show(block=False)
