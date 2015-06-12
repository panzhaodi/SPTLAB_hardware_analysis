import pylab as pl
from glob import glob
import cPickle as pkl
import os, re, bisect, sys
import numpy as np
from scipy.optimize import curve_fit
import sptpol_software.util.time as sptTime
import sptpol_software.util.mpfit as mpfit
from uofc.tools.leastsqbound import leastsqbound
pl.rcParams['legend.fontsize']='small'

# straight up stolen from JT Sayre.
#    - now slightly edited by T.Natoli

Fs = 25e6/2**17
rtdir = '/home/cryo/Data/output/20131220_cooldown/br103/'
imag = np.complex(0,1)

# Below are specific patters we will search log files for
tm_tmpl = '(\d{4}-\d{2}-\d{2})-(\d+_\d+_\d+)'
tmp = '(\d{4})-(\d{2})-(\d{2})-(\d+)_(\d+)_(\d+)'

saved_pngs = []

def get_etf_data_points(data,SQ,ch):
    '''
    this function takes the fft of the raw etf data and finds the amplitude of the 
       24Hz bin as well as an estimate on that amplitude and returns a bunch of stuff, 
       including those two numbers along with the actual fft and raw data
    
    data = data from etf data taken at a particular Rfrac
    SQ = squid
    ch = channel 
    '''
    fq = data['frequencies']
    if not pl.asarray([('Hz' in x) for x in data.keys() if x != 'frequencies']).all():
        ordered_dat = sorted([((int(x.split(' ')[0]) if 'Hz' in x else 100+int(x[1:3])) , x) for x in data.keys() if x!='frequencies'])
    else:
        ordered_dat = sorted([(float(x.split(' ')[0]), x) for x in data.keys() if x!='frequencies'])
    thisDat = {}
    thisErr = {}
    outDict = dict({})

    for i,F in enumerate(fq):
        try:
            dat_i = data[ordered_dat[i][1]]
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

def single_body_model(w,tau):#,norm_freq=None):
    '''
    This model will be normalized to the first given frequency
    '''
    return abs((1+imag*w[0]*tau)/(1+imag*w*tau))

def single_body_model_mp(tau,w,data,err,fjac=None):
    output = (data - single_body_model(w,tau))/err
    status = 0
    return status, output

def body1_lr_model(w, tau_o, Rnow, inductance):
    tau_e = 2.0*inductance/(Rnow)
    model = 1./(np.sqrt(1+w**2*tau_e**2)) * (2*imag*tau_e*tau/tau_o + imag*w*tau - w**2*tau_e*tau + 1 - imag*w*tau_e)**(-1)
    return model


def two_body_model(w,tau,gamma,tau_0):
    model = abs(1.0/( (1+imag*w*tau_0)/(1+imag*w*tau_0/(1+gamma)) + tau_0/tau -1))
    norm_model = abs(model*( (1+imag*w[0]*tau_0)/(1+imag*w[0]*tau_0/(1+gamma)) + tau_0/tau -1))
    return norm_model

def two_body_model_mp(params,w,data,err,fjac=None):
    tau,gamma,tau_0 = params
    output = (data - two_body_model(w,tau,gamma,tau_0))/err
    return 0, output


def two_body_model_lr(w,tau,gamma,tau_0,Rnow,inductance):
    tau_e = 2.0*inductance/(Rnow)
#    model = abs(1.0/( (1+imag*w*tau_0)/((1+imag*w*tau_0/(1+gamma))*(1+imag*w*tau_e)) + (tau_0/tau -1)*(1-imag*w*tau_e)))
 #   norm_model = abs(model*( (1+imag*w[0]*tau_0)/(1+imag*w[0]*tau_0/(1+gamma)) + tau_0/tau -1))##########
    model = lambda w,gamma,tau,tau_0,tau_e: abs((1./np.sqrt(1+w**2*tau_e**2))*(1./(((1+imag*w*tau_0)*(1+imag*w*tau_e))/(1+imag*w*tau_0/(1+gamma))+(tau_0/tau -1)*(1-imag*2*tau_e))))
    norm_model = model(w,gamma,tau,tau_0,tau_e) / model(w[0],gamma,tau,tau_0,tau_e)
    return norm_model

def two_body_model_lr_mp(params,w,data,err,fjac=None):
    tau,gamma,tau_0,Rnow,inductance = params
    output = (data-two_body_model_lr(w,gamma,tau,tau_0,Rnow,inductance))/err
    return 0, output
    
def fit_data(fit_type,freqs,amps,errs,
               f_3db=10.0,
               gamma = 300,
               tau_0 = 60E-3,
               Rnow = 0,
               inductance = 0):
    ''' fit the specified model

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
    if fit_type == 'single_body':
        p0 = [ tau ]# guess at initial tau
        functkw = {'w':w,'data':amps, 'err':errs}
        fit_output = mpfit.mpfit(single_body_model_mp, p0, 
                                 functkw = functkw,quiet=True)
        tau,tau_err = fit_output.params[0],fit_output.perror[0]
        print asdf
        return tau,tau_err

    elif fit_type == 'two_body':
        p0 = [tau,gamma,tau_0]
        functkw = {'w':w,'data':amps,'err':errs}
        fit_output = mpfit.mpfit(two_body_model_mp, p0, 
                                 functkw = functkw,quiet=False)
        tau,gamma,tau_0 = fit_output.params
        tau_err,gamma_err,tau_0_err = fit_output.perror
        
        return abs(tau), tau_err, gamma, gamma_err, tau_0, tau_0_err

    elif fit_type == 'two_body_lr':
        if not inductance:
            inductance = 20E-6 # Henries
            print 'You did not supply an INDUCTANCE, I am using inducatance ='+str(inductance)+' Henries'
        if not Rnow:
            Rnow = 0.261 #Ohms
            print 'You did not supply a RESISTANCE to the fitter, I am using Rnow ='+str(Rnow)+' Ohms'
        
        p0 = [tau,gamma,tau_0,Rnow,inductance]
        functkw = {'w':w,'data':amps,'err':errs}
        # freeze the inductance and resistance
        parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for i in range(len(p0))] 
        for i in range(len(p0)): parinfo[i]['value']=p0[i] 
#        parinfo[1]['limited'] = [1,1]
#        parinfo[1]['limits'] = [0,2000] # limits the values of gamma
        parinfo[1]['fixed'] = 1 # freezes the gamma
        parinfo[3]['fixed'] = 1 # freezes the Rnow
        parinfo[4]['fixed'] = 1 # freezes the inductance

        fit_output = mpfit.mpfit(two_body_model_lr_mp, p0, parinfo=parinfo,
                                 functkw = functkw,quiet=False)
        tau,gamma,tau_0,Rnow,inductance =  fit_output.params
        print fit_output.params
        tau_err,gamma_err,tau_0_err,Rnow_err,inductance_err =  fit_output.perror
        return abs(tau), tau_err, gamma, gamma_err, tau_0, tau_0_err
    else:
        print 'You did not supply a correct model option to the fitter'
        return abs(tau), tau_err, gamma, gamma_err, tau_0, tau_0_err

def make_plots(data,savePlots=True,directory='',
               single_body = False, two_body = False, 
               two_body_lr = False,
               fit_freq_cutoff = 1000,
               inductance = 0,
               Rt = 0,
               gamma = 200,
               verbose=False):
    '''
    In addition to making visually pleasing plots, 
       this function calls the fitting routine for the data.

    '''
    for i,bolo in enumerate(data.keys()):
        _now = data[bolo]
        #initiate a figure to use
        pl.figure(i+1,figsize=(12,10));pl.clf();
        pl.plot(1,1,'k',linestyle='none',label= '   %s   %s  %s %s'%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(ms)$',r'$\gamma$'))
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


            # plot the data
            pl.errorbar(freqs,respnorm,yerr=errnorm,fmt = '.',#label=label,
                        color = pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))

            ###################
            # we probably only want to use frequencies aboave a certain point 
            print 'Only fitting frequencies below %s.'%(str(fit_freq_cutoff))
            ind_to_fit = pl.find(freqs < fit_freq_cutoff)[-1]
#                low_ind_to_fit = pl.find(freqs >= fit_freq_min)[0]
#            ind_to_fit = -1
            # only plot the fit line to the frequencies we fit to
            #freqs_plot = pl.linspace(freqs[0],freqs[ind_to_fit],500)# a dummy list to plot a smooth fit line
            freqs_plot = pl.logspace(1,np.log(freqs[-1])/np.log(freqs[0]),base = freqs[0])
#            pl.loglog(freqs_plot,single_pole(2*pl.pi*freqs_plot,tau),color=pl.gca().lines[-1].get_color(),linewidth=0.4)
            if single_body:
                tau,tau_err = fit_data('single_body',freqs[:ind_to_fit],
                                        respnorm[:ind_to_fit],
                                        errnorm[:ind_to_fit], f_3db = f_3db)
                gamma = 'inf'
                model_to_plot = single_body_model(2*pl.pi*freqs_plot,tau)

            if two_body:
                fit_output = fit_data('two_body',freqs[:ind_to_fit],
                                        respnorm[:ind_to_fit],
                                        errnorm[:ind_to_fit], f_3db = f_3db)
                tau,tau_err,gamma,gamma_err,tau_0,tau_0_err = fit_output
                model_to_plot = two_body_model(2*pl.pi*freqs_plot,tau,gamma,tau_0)

            if two_body_lr:
                fit_output = fit_data('two_body_lr',freqs[:ind_to_fit],
                                      respnorm[:ind_to_fit],
                                      errnorm[:ind_to_fit], f_3db = f_3db,
                                      inductance = inductance, Rnow = Rt*float(rfrac),
                                      gamma = gamma)
                tau,tau_err,gamma,gamma_err,tau_0,tau_0_err = fit_output
                model_to_plot = two_body_model_lr(2*pl.pi*freqs_plot,tau,gamma,tau_0,Rt*float(rfrac),inductance)



            # if tau is unresponably high (corresponds to a timeconstant >500)
            #    set tau = NA in the label
            if 1/(2*pl.pi*tau) < 500:
                label = '%8.3f  %6.1f  %6.1f(%6.2f)  %6.1f(%6.2f)'%(float(rfrac),1/(2*pl.pi*tau),tau*1000.,tau_err*1000.,gamma,gamma_err)
            else:
                label = '%8.3f  %s  %6.1f'%(float(rfrac),'  NA  ',tau*1000.)
            print label,'+/-', tau_err*1000

            pl.loglog(freqs_plot,model_to_plot,color=pl.gca().lines[-1].get_color(),linewidth=0.4,label = label)

        pl.yscale('log'); pl.xscale('log')
        pl.legend()
        pl.xlabel('freq (Hz)')
        pl.ylabel('Response (normalized to lowest freq)')
        pl.suptitle('%s'%directory)
        pl.title('Quick Look ETF (%s)'%(bolo))
        pl.xlim(freqs[0],freqs[-1]);pl.ylim(1e-3,2)
        pl.grid()
        if savePlots:
            pl.savefig(rtdir + '%s/TauETF_MPIT_err_%s'%(directory,bolo))
            saved_pngs.append(rtdir+'%s/TauETF_MPIT_err_%s.png'%(directory,bolo))



  
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

    try:
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

    for _sub in subs:
        analyze_etf_dir(_sub)

def analyze_etf_dir(directory, save=False, Rn = None,
                    single_body=True,
                    two_body=False,
                    two_body_lr=False,
                    fit_freq_cutoff = 1000,
                    inductance = 20E-6,
                    gamma = 200,
                    Rt = 0.261):
    '''  analyzes all of the etf data in a given directory
    '''
    
    print 'starting analysis of ' + directory
    data ={}
    # get the squid names from the pkl file names
    all_etf_pkls = glob('%s/*_TauETF*.pkl'%(directory))
    squids = [_pkl.rpartition('/')[-1].partition('_')[0] for _pkl in all_etf_pkls]
    squids = list(set(squids))
    totalDict=dict()

    for _SQ in squids:
            
        # gather and extract the timestamp from the tauETF pkl and bolos_tune logs
        tau_pkls = glob('%s/%s_TauETF*.pkl'%(directory,_SQ))
        tau_times =[sptTime.toSptDatetime(_pkl.rpartition('_')[-1].rpartition('.pkl')[0]).fil 
                    for _pkl in tau_pkls]
        
#        tune_pkls = glob('%s/%s*_TuneBoloCombDan_*.pkl'%(directory,_SQ))

        tune_logs = glob('%s/AlgOutTuneBoloCombDan*.log'%(directory))
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
            print _pkl
            tunes_index = bisect.bisect(tune_times,_time)-1
            # reach into the bolo_tune.log to find the Rfrac
            rFrac,success = get_Rfrac_from_tune_log(tune_logs[tunes_index],_SQ,directory)
            if not success:
                break
            # open the etf pickle and fit 
            thisFL = pkl.load(open(_pkl))
            for ch in [x for x in thisFL.keys() if 'Ch' in x]:
                data_points,err,outDict = get_etf_data_points(thisFL[ch],_SQ,ch)
                try:
                    complete_crap = totalDict[_SQ+ch]
                except KeyError:
                    totalDict[_SQ+ch] = dict({rFrac[ch]:dict()})
                totalDict[_SQ+ch][rFrac[ch]]={'tickle_frequencies':outDict}
                data.setdefault(_SQ+ch,{})
                data[_SQ+ch][rFrac[ch]] = data_points

        make_plots(totalDict,savePlots=save,directory=directory,
                   single_body=single_body,
                   two_body=two_body,
                   two_body_lr=two_body_lr,
                   fit_freq_cutoff=fit_freq_cutoff,
                   Rt=Rt,
                   gamma = gamma)
        if save:
            pkl.dump(totalDict,open( directory + '/TauETF_all.pkl','w'))
    
    print 'Sucessfully saved %i pngs' %(len(saved_pngs))
    for _png in saved_pngs:
        print _png
                                










##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
        



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
