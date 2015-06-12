import pylab as pl
from glob import glob
import cPickle as pkl
import os, re, bisect, sys
import numpy as np
from scipy.optimize import curve_fit
import sptpol_software.util.time as sptTime
pl.rcParams['legend.fontsize']='small'

# straight up stolen from JT Sayre.

Fs = 25e6/2**17
rtdir = '/home/cryo/Data/output/20140414_cooldown/'
imag = np.complex(0,1)
#subs = [x for x in os.listdir(rtdir) if os.path.isdir(rtdir+x)]
#subs = ['20140111_224853_IV_ETF']
#squids = ['Sq8']#,'Sq8']

tm_tmpl = '(\d{4}-\d{2}-\d{2})-(\d+_\d+_\d+)'
tmp = '(\d{4})-(\d{2})-(\d{2})-(\d+)_(\d+)_(\d+)'
saved_pngs = []

def calc_dumb_etf(s1,sb,SQ,ch,rFrac,test=False):
    fq = s1['frequencies']
    if not pl.asarray([('Hz' in x) for x in s1.keys() if x != 'frequencies']).all():
        ordered_dat = sorted([((int(x.split(' ')[0]) if 'Hz' in x else 100+int(x[1:3])) , x) for x in s1.keys() if x!='frequencies'])
    else:
        ordered_dat = sorted([(float(x.split(' ')[0]), x) for x in s1.keys() if x!='frequencies'])
    thisDat = {}
    thisErr = {}
    outDict = dict({})
#    pl.figure(20,figsize=(12,10));pl.clf();
    for i,F in enumerate(fq):
        try:
            dat_i = s1[ordered_dat[i][1]]
        except IndexError:
            print 'The frequencies are mismatched for '+_pkl
            print '***I will continue, but the plot will not make any sense!!'
            break
        N = len(dat_i)
        hanning_window = pl.hanning(N)
        FT, FQ = abs(pl.fft(dat_i*hanning_window))[:N/2], pl.fftfreq(N,1/Fs)[:N/2]
        pkind = pl.argmin(abs(FQ-24))
#        max_amp = FT[pkind]
#        pl.find(FT>0.1*max_amp)
        pkl_quarter_inds = pl.find(abs(FQ-24)<0.25)

        if len(pkl_quarter_inds)>3:
            thisDat['%.1f'%F] = sum(FT[pkl_quarter_inds])
            thisErr['%.1f'%F] = (sum(FT[pkl_quarter_inds+len(pkl_quarter_inds)])+
                                 sum(FT[pkl_quarter_inds-len(pkl_quarter_inds)]))/2
        else:
            thisDat['%.1f'%F] = sum(FT[(pkind-1):(pkind+2)])
            thisErr['%.1f'%F] = (sum(FT[pkind+3:pkind+6])+
                                 sum(FT[pkind-5:pkind-2]))/2

#        if i >4 :
#            print asdf
        outDict[F] = dict({'fft':FT,'fft_freq':FQ-24,
                           'TOD':dat_i,'final_value':thisDat['%.1f'%F],
                           'value_error':thisErr['%.1f'%F],
                           })
        
    return thisDat,thisErr,outDict

def single_pole(w,tau):#,norm_freq=None):
#    return 1/pl.sqrt(1+(w*tau)**2)
#    return 1/abs(1+imag*w*tau)
    '''
    norm_freq = Frequency to normalize the single pole at. 
               If None, then the first frequency given will be used.
    '''
#    if norm_freq:
#        return abs((1+imag*norm_freq*tau)/(1+imag*w*tau))
#    else:
#        print 'this anyways'
    return abs((1+imag*w[0]*tau)/(1+imag*w*tau))
    

def body1_lr_model(w, tau_o, Rn, Rfrac, inductance = 20):
    tau_e = 2.0*inductance/(Rfrac*Rn)
    
    model = 1./(np.sqrt(1+w**2*tau_e**2)) * (2*imag*tau_e*tau/tau_o + imag*w*tau - w**2*tau_e*tau 
                                             + 1 - imag*w*tau_e)**(-1)
    return model

def two_pole_model(w,gamma,tau,tau_0):
#    tau_e = 2.0*inductance/(Rfrac*Rn)
    model = abs(1.0/( (1+imag*w*tau_0)/(1+imag*w*tau_0/(1+gamma)) + tau_0/tau -1))
    norm_model = abs(model*( (1+imag*w[0]*tau_0)/(1+imag*w[0]*tau_0/(1+gamma)) + tau_0/tau -1))

    return norm_model

def two_pole_model_lr(w,gamma,tau,tau_0,Rt):
    inductance = 20E-6# Henries
    tau_e = 2.0*inductance/(Rt)
    model = abs(1.0/( (1+imag*w*tau_0)/((1+imag*w*tau_0/(1+gamma))(1+imag*w*tau_e)) + (tau_0/tau -1)*(1-imag*w*tau_e)))
    norm_model = abs(model*( (1+imag*w[0]*tau_0)/(1+imag*w[0]*tau_0/(1+gamma)) + tau_0/tau -1))##########

    return norm_model


def single_pole_fit(freqs,amps,errs,f_3db=10.0):
    ''' just fit a single pole to the data given

    inputs: 
        freqs: a list of frequencies 
        amps: a list of amplitudes corresponding to the freqs
        tau_0: a guess at the tau of the pole you are fitting
    outputs:
        the fit tau
    '''
    
    w = 2*pl.pi*pl.asarray(freqs) # convert frequency into angular frequency 
    tau = float(f_3db)/(2*pl.pi*f_3db)
#    popt,pcov = curve_fit(single_pole,w,amps,tau_0)
    popt,pcov = curve_fit(single_pole,w,amps,tau,sigma=errs)
#    print asdf
    if not np.isfinite(pcov):
        perr = np.inf
    else:
        perr = np.sqrt(pcov[0][0])
    return abs(popt[0]),perr

def make_plots(data,sb,savePlots=True):

    for sb in data.keys():

        
        for i,bolo in enumerate(data[sb].keys()):
            _now = data[sb][bolo]
            #initiate a figure to use
            pl.figure(i+1,figsize=(12,10));pl.clf();
            pl.plot(1,1,'k',linestyle='none',label= '   %s   %s  %s  %s'%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(ms)$',r'$\pm\,\tau\~(ms)$'))
            # cycle through the responses at different Rfractions
            rfrac_count = len(_now.keys())
            for _nrfrac, rfrac in enumerate(sorted(_now.keys(),reverse=True)):
                freqs,responses,errors = [],[],[]
                for _freq in _now[rfrac]['tickle_frequencies'].keys():
                    _now_tick = _now[rfrac]['tickle_frequencies'][_freq]
                    freqs.append(float(_freq))
                    responses.append(_now_tick['final_value'])
                    errors.append(_now_tick['value_error'])
                    # reorder the frequencies and responses to increasing in freq
                freqs,responses,errors = pl.array(zip(*sorted(zip(freqs,responses,errors))))
                
                # normalize the response to the first point
                respnorm = pl.asarray(responses)/responses[0]
                errnorm = pl.asarray(errors)/responses[0]
                # finding the 3db point? T-funk, check this out more
                ii = min([pl.find(respnorm>0.5)[-1], len(respnorm)-2])    ### check this
                f_3db = pl.interp(pl.sqrt(0.5),respnorm[[ii+1,ii]],freqs[[ii+1,ii]])  ### check this
                pl.axhline(pl.sqrt(0.5),color='k',linestyle='--')
            ###################
            # we probably only want to use frequencies aboave a certain point 
            #   right now we are hard coding it to only fit with freqs below 1,000Hz
                fit_freq_cutoff = 1000.0 # in Hz
                #fit_freq_min = 1.0
                print 'Only fitting frequencies below %s.'%(str(fit_freq_cutoff))
                ind_to_fit = pl.find(freqs < fit_freq_cutoff)[-1]
#                low_ind_to_fit = pl.find(freqs >= fit_freq_min)[0]
#            ind_to_fit = -1
                tau,tau_err = single_pole_fit(freqs[:ind_to_fit],
                                              respnorm[:ind_to_fit],
                                              errnorm[:ind_to_fit], f_3db = f_3db)

            # if tau is unresponably high (corresponds to a timeconstant >500)
            #    set tau = NA in the label
                if 1/(2*pl.pi*tau) < 500:
                    label = '%8.3f  %6.1f  %6.1f  %6.2f'%(float(rfrac),1/(2*pl.pi*tau),tau*1000.,tau_err*1000.)
                else:
                    label = '%8.3f  %s  %6.1f'%(float(rfrac),'  NA  ',tau*1000.)
                print label,'+/-', tau_err*1000

#                pl.loglog(freqs,respnorm,'.',markersize=12,label=label,color = pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))
                pl.errorbar(freqs,respnorm,yerr=errnorm,fmt = '.',label=label,
                            color = pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))
#                            markeredgecolor='k')
#
            # only plot the fit line to the frequencies we fit to
                freqs_plot = pl.linspace(freqs[0],freqs[ind_to_fit],500)# a dummy list to plot a smooth fit line
            #freqs_plot = pl.linspace(freqs[0],freqs[-1],500)# To plot fit across all freqs
                pl.loglog(freqs_plot,single_pole(2*pl.pi*freqs_plot,tau),color=pl.gca().lines[-1].get_color(),linewidth=0.4)
            pl.yscale('log'); pl.xscale('log')
            pl.legend()
            pl.xlabel('freq (Hz)')
            pl.ylabel('Response (normalized to lowest freq)')
            pl.suptitle('%s'%sb)
            pl.title('Quick Look ETF (%s)'%(bolo))
            pl.xlim(freqs[0],freqs[-1]);pl.ylim(1e-3,2)
            pl.grid()
            if savePlots:
                pl.savefig(rtdir + '%s/TauETF_MPIT_err_%s'%(sb,bolo))
                saved_pngs.append(rtdir+'%s/TauETF_MPIT_err_%s.png'%(sb,bolo))



  
def get_Rfrac_from_tune_log(thisLog,SQ,sb):
    # take a look inside the bolo_tune.log to find the Rfrac
    #
    #This assumes all bolos were biased at the same Rn, so the log will list the same number for each channel
    tmpl = "%s: \{('channel_frac_resistance': \[([\d\.\s,]+)+\]).*'channels': \[([\d,\s]+)+\]"%SQ
    out={}
    success = True
    with open(thisLog) as log_file:
        for ln in log_file:
            # find the line we want, the one that starts the IV
            #  and reports the resistance
            if SQ in ln and 'channel_frac_resistance' in ln:
                resvals, chans = pl.asarray(re.search(tmpl,ln).groups())[1:3]
                resvals,chans = resvals.split(','), chans.split(',')
            if SQ in ln and 'Traceback' in ln:
                print 'Error in ' + thisLog + ' IV curve, not analyzing'
                success = False
                
    return dict([('Ch'+x.strip(),y.strip()) for x,y in zip(chans,resvals)]), success

if __name__ == '__main__':

    test = False
    try:
        if 'test' in sys.argv:
            print asdf
            test = True
        else:
            subs = sys.argv[1:]
        print subs
    except  IndexError:
        print "  How to run:"
        print "   >>python plot_etf.py directory"
        sys.exit()
    else:
        pass

    if sys.argv[1] == "help":
        print "   How to run:"
        print "   >>python plot_etf.py directory "
        sys.exit()
    
    totalDict=dict()

    for sb in subs:
        print sb
        data ={}
        # get the squid names from the pkl files instead of listing them here
        all_etf_pkls = glob(rtdir+ '%s/*_TauETF*.pkl'%(sb))
        squids = [_pkl.rpartition('/')[-1].partition('_')[0] for _pkl in all_etf_pkls]
        squids = list(set(squids))
        totalDict[sb]=dict()
#        for SQ in ['SB1Ch2','SB1Ch4','SB1Ch6','SB1Ch8']:
        for SQ in squids:
            
            ky = SQ
            # gather and extract the timestamp from the tauETF pkl and bolos_tune logs
            tau_pkls = glob(rtdir+ '%s/%s_TauETF*.pkl'%(sb,SQ))
            tau_times =[sptTime.toSptDatetime(_pkl.rpartition('_')[-1].rpartition('.pkl')[0]).fil 
                        for _pkl in tau_pkls]
            tune_logs = glob(os.path.join(rtdir,sb)+'/AlgOutTuneBoloCombDan*.log')
            tune_times = [sptTime.toSptDatetime(_log.rpartition('Dan')[-1].rpartition('.log')[0]).fil 
                          for _log in tune_logs]
            
            # sort the lists we just got by time
            taus = zip(tau_times, tau_pkls)
            taus.sort()
            tau_times, tau_pkls = zip(*taus)
            tunes = zip(tune_times, tune_logs)
            tunes.sort()
            tune_times, tune_logs = zip(*tunes)

            for (_time,_pkl) in taus:
                print _time, _pkl
                tunes_index = bisect.bisect(tune_times,_time)-1
                # reach into the bolo_tune.log to find the Rfrac
                rFrac,success = get_Rfrac_from_tune_log(tune_logs[tunes_index],SQ,sb)
                if not success:
                    break

                # open the etf pickle and fit 
                thisFL = pkl.load(open(_pkl))
                for ch in [x for x in thisFL.keys() if 'Ch' in x]:
                    dat,err,outDict = calc_dumb_etf(thisFL[ch],sb,SQ,ch,rFrac,test=test)
                    try:
                        complete_crap = totalDict[sb][SQ+ch]
                    except KeyError:
                        totalDict[sb][SQ+ch] = dict({rFrac[ch]:dict()})
                    totalDict[sb][SQ+ch][rFrac[ch]]={'tickle_frequencies':outDict}
                    data.setdefault(ky+ch,{})
                    data[ky+ch][rFrac[ch]] = dat

#        make_plots_old(data,sb)
        make_plots(totalDict,sb)
        pkl.dump(totalDict,open(rtdir +'/'+ sb + '/TauETF_all.pkl','w'))
    
    print 'Sucessfully saved %i pngs' %(len(saved_pngs))
    for _png in saved_pngs:
        print _png
                                













def investigate(data, savePlots=False,SqCh=None,Rfrac=None,
                fit_freq_cutoff = 1000,extras = [],checkTau = False,
                checkGamma = False, tau_0 = 50E-3):

    for sb in data.keys():
        if SqCh:
            SqChs = [SqCh]
        else:
            squids = data[sb].keys()
        for i,bolo in enumerate(SqChs):
            _now = data[sb][bolo]
            #initiate a figure to use
            pl.figure(i+1,figsize=(12,10));pl.clf();
            axes = pl.subplot(111)
            pl.plot(1,1,'k',linestyle='none',label= '   %s   %s  %s  %s '%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(ms)$',r'$\pm\~(ms)$'))
#            main_label_info = '   %s   %s  %s  %s '%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(ms)$',r'$\pm\~(ms)$')

            # cycle through the responses at different Rfractions
            if Rfrac:
                try:
                    total_crap = _now[Rfrac]
                except KeyError:
                    try:
                        total_crap = _now[str(Rfrac)]
                        Rfrac = str(Rfrac)
                    except KeyError:
                        print 'You entered an Rfrac that is not in the data.'
                        print 'Availible Rfracs are :'+str(_now.keys())
                        return
                rfracs = [Rfrac]
            else:
                rfracs = _now.keys()
            rfrac_count = len(rfracs)
            for _nrfrac, rfrac in enumerate(sorted(rfracs,reverse=True)):
                freqs,responses,errors = [],[],[]
                for _freq in _now[rfrac]['tickle_frequencies'].keys():
                    _now_tick = _now[rfrac]['tickle_frequencies'][_freq]
                    freqs.append(float(_freq))
                    responses.append(_now_tick['final_value'])
                    errors.append(_now_tick['value_error'])
                    # reorder the frequencies and responses to increasing in freq
                freqs,responses,errors = pl.array(zip(*sorted(zip(freqs,responses,errors))))
                # normalize the response to the first point
                respnorm = pl.asarray(responses)/responses[0]
                errnorm = pl.asarray(errors)/responses[0]
                # finding the 3db point? T-funk, check this out more
                ii = min([pl.find(respnorm>0.5)[-1], len(respnorm)-2])    ### check this
                f_3db = pl.interp(pl.sqrt(0.5),respnorm[[ii+1,ii]],freqs[[ii+1,ii]])  ### check this
                pl.axhline(pl.sqrt(0.5),color='k',linestyle='--')
            ###################
            # we probably only want to use frequencies aboave a certain point 
            #   right now we are hard coding it to only fit with freqs below fit_freq_cutoff
                print 'Only fitting frequencies below %s'%(str(fit_freq_cutoff))
                ind_to_fit = pl.find(freqs < fit_freq_cutoff)[-1]
#            ind_to_fit = -1
                tau, tau_err = single_pole_fit(freqs[:ind_to_fit],respnorm[:ind_to_fit],
                                               errnorm[:ind_to_fit],f_3db)
                
            # if tau is unresponably high (corresponds to a timeconstant >500)
            #    set tau = NA in the label
                if 1/(2*pl.pi*tau) < 500:
                    main_label = '%8.3f  %6.1f  %6.1f  %6.2f'%(float(rfrac),1/(2*pl.pi*tau),tau*1000.,tau_err*1000)
                else:
                    main_label = '%8.3f  %s  %6.1f'%(float(rfrac),'    NA   ',tau*1000.)
                print main_label, '+/-', tau_err*1000

#                pl.loglog(freqs,respnorm,'.',markersize=12,label=label,color = pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))
                pl.errorbar(freqs,respnorm,yerr=errnorm,fmt = '.',label=main_label,color = 'b')#pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))
                
            # only plot the fit line to the frequencies we fit to
#                freqs_plot = pl.linspace(freqs[0],freqs[ind_to_fit],1000)# a dummy list to plot a smooth fit line
#                freqs_plot = pl.linspace(freqs[0],freqs[-1],1000)# To plot fit across all freqs
                freqs_plot = [freqs[0]]
                while freqs_plot[-1] < 10000:
                    freqs_plot.append(freqs_plot[-1]*1.1)
                freqs_plot = np.array(freqs_plot)
                pl.loglog(freqs_plot,single_pole(2*pl.pi*freqs_plot,tau),color=pl.gca().lines[-1].get_color(),linewidth=0.8)

                colors = ['g','m','k','y','b','g','m','k','g','y','b','g','m','k','g','y']
                colors = colors[:len(extras)]
                if checkTau:
                    pl.plot(1,1,'k',linestyle='none',label= '   %s   %s  %s  %s '%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(ms)$',r'$\pm\~(ms)$'))
                    for _etau,_color in zip(extras,colors):
                        pl.loglog(freqs_plot,single_pole(2*pl.pi*freqs_plot,_etau),linewidth=0.5,color=_color,
                                  label = '  %s  %6.1f  %6.1f'%('  NA  ',1/(2*pl.pi*_etau),_etau*1000.))
                elif checkGamma:
                    pl.plot(1,1,'k',linestyle='none',label= '   %s   %s  %s  %s '%(r' $\gamma$ ',r' $\tau_{0}\,(ms)$ ',r'  $\tau\,(ms)$ ',r'$\pm\~(ms)$'))
                    for _gamma,_color in zip(extras,colors):
                        pl.loglog(freqs_plot,two_pole_model(2*pl.pi*freqs_plot,_gamma,tau,tau_0),linewidth=0.5,color=_color,
                                  label = '%6.0f   %6.1f    %6.1f'%(_gamma,tau_0*1000.,tau*1000.))
            pl.yscale('log'); pl.xscale('log')
            legend_handles,legend_labels = axes.get_legend_handles_labels()

            new_order = [0,-1]
            new_order.extend(range(1,len(legend_labels)-1))
            new_order = np.array(new_order)
#            pl.legend(numpoints = 1,loc=1)
            pl.legend(np.array(legend_handles)[new_order],np.array(legend_labels)[new_order],numpoints = 1,loc=1)

            pl.xlabel('freq (Hz)')
            pl.ylabel('Response (normalized to lowest freq)')
            pl.suptitle('%s'%(sb))
            pl.title('Investigate ETF (%s)'%(bolo))
            pl.xlim(freqs[0],freqs[-1]);pl.ylim(1e-3,2)
            pl.grid()
            if savePlots:
                pl.savefig(rtdir + '%s/TauETF_MPIT_err_%s'%(sb,bolo))
                saved_pngs.append(rtdir+'%s/TauETF_MPIT_err_%s.png'%(sb,bolo))























































def make_plots_old(data,sb,inductance=20,savePlots=False):
    '''
    This is old, do not use, it is just for posterity now
    '''

    for i,bolo in enumerate(data.keys()):
        #initiate a figure to use
        pl.figure(i+1,figsize=(12,10));pl.clf();
        pl.plot(1,1,'k',linestyle='none',label= '   %s   %s  %s'%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(ms)$'))
        # cycle through the responses at different Rfractions
        rfrac_count = len(data[bolo].keys())
        for _nrfrac, rfrac in enumerate(sorted(data[bolo].keys(),reverse=True)):
            freqs,responses = [],[]
            for _freq,_response in data[bolo][rfrac].iteritems():
                freqs.append(float(_freq))
                responses.append(_response)
            # reorder the frequencies and responses to increasing in freq
            freqs,resp = pl.array(zip(*sorted(zip(freqs,responses))))
            # normalize the response to the first point
            respnorm = pl.asarray(resp)/resp[0]
            # finding the 3db point? T-funk, check this out more
            ii = min([pl.find(respnorm>0.5)[-1], len(respnorm)-2])    ### check this
            f_3db = pl.interp(pl.sqrt(0.5),respnorm[[ii+1,ii]],freqs[[ii+1,ii]])  ### check this
            pl.axhline(pl.sqrt(0.5),color='k',linestyle='--')
            ###################
            # we probably only want to use frequencies aboave a certain point 
            #   right now we are hard coding it to only fit with freqs below 1,000Hz
            fit_freq_cutoff = 1000 # in Hz
            print 'Only fitting frequencies below %s'%(str(fit_freq_cutoff))
            ind_to_fit = pl.find(freqs < fit_freq_cutoff)[-1]
#            ind_to_fit = -1
            tau = single_pole_fit(freqs[:ind_to_fit],respnorm[:ind_to_fit],f_3db)

            # if tau is unresponably high (corresponds to a timeconstant >500)
            #    set tau = NA in the label
            if 1/(2*pl.pi*tau) < 500:
                label = '%8.3f  %6.1f  %6.1f'%(float(rfrac),1/(2*pl.pi*tau),tau*1000.)
            else:
                label = '%8.3f  %s  %6.1f'%(float(rfrac),'  NA  ',tau*1000.)

            pl.loglog(freqs,respnorm,'.',markersize=12,label=label,color = pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))
            # only plot the fit line to the frequencies we fit to
            freqs_plot = pl.linspace(freqs[0],freqs[ind_to_fit],500)# a dummy list to plot a smooth fit line
            #freqs_plot = pl.linspace(freqs[0],freqs[-1],500)# To plot fit across all freqs
            pl.loglog(freqs_plot,single_pole(2*pl.pi*freqs_plot,tau),color=pl.gca().lines[-1].get_color(),linewidth=0.4)
        pl.legend()
        pl.xlabel('freq (Hz)')
        pl.ylabel('Response (normalized to lowest freq)')
        pl.title('Quick Look ETF (%s)'%bolo)
        pl.xlim(freqs[0],freqs[-1]);pl.ylim(1e-3,1e1)
        pl.grid()
#        pl.show()
        # savefig('/home/sayre/spider/3g_testing/20131222_025713/%s_TauETF_MPIT'%bolo)
        if savePlots:
            pl.savefig(rtdir + '%s/TauETF_MPIT_%s'%(sb,bolo))
            saved_pngs.append(rtdir+'%s/TauETF_MPIT_%s.png'%(sb,bolo))
