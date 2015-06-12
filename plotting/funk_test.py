import pylab as pl
import numpy as np






'''

etf_data = zip([0.995,0.990,0.960,0.930,0.900,0.870,0.800,0.700],
               [35.2 ,16.0 ,9.5  ,6.9  ,6.8  ,6.6  ,5.7  ,14.0])
iv_data = [ivs_p1,ivs_p09,ivs_p08]


'''


def tau_0_plots(iv_data,etf_data):

    _sq = 'Sq8Ch1'

    pl.figure(30)
    pl.clf()
    for _RparDict in iv_data:
        
        _Rpar = _RparDict[_sq]['Rparasitic']
        Rfracs = np.sort(_RparDict[_sq]['Rfrac'].keys())[::-1]
        etf_Rfracs,taus_raw = zip(*etf_data)
        Rfracs = list(set(Rfracs).intersection(etf_Rfracs))
        loopgains = np.array([_RparDict[_sq]['Rfrac'][_rf]['loopgain'][0] for _rf in Rfracs])
        final_smooth_loopgains = np.array([_RparDict[_sq]['Rfrac'][_rf]['final_smooth_loopgain'] for _rf in Rfracs])

        
        taus = np.array([taus_raw[etf_Rfracs.index(_rf)] for _rf in Rfracs])
        pl.plot(Rfracs,taus*(loopgains+1),'o',label=str(_Rpar))
        pl.plot(Rfracs,taus*(final_smooth_loopgains+1),'^',label=str(_Rpar)+' smooth',
                color = pl.gca().lines[-1].get_color())

    pl.legend(title = 'Rpar')
    pl.grid()
    pl.ylabel('Tau_0 [ms] (= tau_etf*(loopgain_iv+1))')
    pl.xlabel('R_frac')
