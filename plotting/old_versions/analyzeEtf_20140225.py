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
#rtdir = '/home/cryo/Data/output/20131220_cooldown/br103/'
imag = np.complex(0,1)

# Below are specific patters we will search log files for
tm_tmpl = '(\d{4}-\d{2}-\d{2})-(\d+_\d+_\d+)'
tmp = '(\d{4})-(\d{2})-(\d{2})-(\d+)_(\d+)_(\d+)'

saved_pngs = []

def get_etf_data_points(data,SQ,ch,plot=False):
    '''
    this function takes the fft of the raw etf data and finds the amplitude of the 
       24Hz bin as well as an estimate on that amplitude and returns a bunch of stuff, 
       including those two numbers along with the actual fft and raw data
    
    data = data from etf data taken at a particular Rfrac
    SQ = squid
    ch = channel 
    '''
    tickle_freq = data['frequencies']
    if not pl.asarray([('Hz' in x) for x in data.keys() if x != 'frequencies']).all():
        ordered_dat = sorted([((int(x.split(' ')[0]) if 'Hz' in x else 100+int(x[1:3])) , x) for x in data.keys() if x!='frequencies'])
    else:
        ordered_dat = sorted([(float(x.split(' ')[0]), x) for x in data.keys() if x!='frequencies'])
    thisDat,thisMeas,thisErr = {},{},{}
    outDict = dict({})
    if plot:
        pl.figure(20)
        pl.clf()

    for i,_tfreq in enumerate(tickle_freq):
        try:
            dat_i = data[ordered_dat[i][1]]
        except IndexError:
            print 'The frequencies are mismatched for '+_pkl
            print '***I will continue, but the plot will not make any sense!!'
            break

        n_pts_in_dat_i = len(dat_i)
        hanning_window = pl.hanning(n_pts_in_dat_i)
        FQ = pl.fftfreq(n_pts_in_dat_i,1/Fs)[:n_pts_in_dat_i/2]
        FT = abs(pl.fft(dat_i*hanning_window))[:n_pts_in_dat_i/2]
        pkind = pl.argmin(abs(FQ-24))

#        max_amp = FT[pkind]
#        pl.find(FT>0.1*max_amp)

        pkl_quarter_inds = pl.find(abs(FQ-24)<0.25)
        if plot:
            pl.plot(FQ-24,FT,'.-',label=str(_tfreq),alpha=0.4)
            pl.plot(FQ[pkl_quarter_inds]-24,FT[pkl_quarter_inds],'o',color = pl.gca().lines[-1].get_color())

        if len(pkl_quarter_inds)>3:
            bins_above = sum(FT[pkl_quarter_inds+len(pkl_quarter_inds)])
            bins_below = sum(FT[pkl_quarter_inds-len(pkl_quarter_inds)])

            now_sig = sum(FT[pkl_quarter_inds])
        else:
            bins_above = sum(FT[pkind+3:pkind+6])
            bins_below = sum(FT[pkind-5:pkind-2])
            now_sig = sum(FT[(pkind-1):(pkind+2)])


        valid = True
        # try:
        #     p0 = [max(FT[pkl_quarter_inds]), 0, 0.2]#Amp, mu, sigma
        #     print asdf
        #     coeff, var_matrix = curve_fit(gauss, FQ[pkl_quarter_inds]-24, FT[pkl_quarter_inds], p0=p0)

        #### ***now_sig = area under the curve
        #     now_sig = coeff[0]*coeff[2]*np.sqrt(2*np.pi)
        # except RuntimeError:
        #     now_sig = sum(FT[(pkind-1):(pkind+2)])
        #     valid = False

        # take out the noise bias from the signal
        if bins_above > now_sig or bins_below > now_sig:
            valid = False
        now_err = (bins_above+bins_below)/2.0
        ######## multiplying the signal by the noise bias correction factor for 
        #####       the S/N measured
#        thisDat['%.1f'%_tfreq] = now_sig*get_bias_correction_factor(now_sig/
#        thisDat['%.1f'%_tfreq] = np.sqrt(now_sig**2-now_err**2)######### I don't think this is the correct thing to do
        print '*****************Not taking out a noise bias**************************'
        thisDat['%.1f'%_tfreq] = now_sig
        thisMeas['%.1f'%_tfreq] = now_sig
        thisErr['%.1f'%_tfreq] = now_err

        outDict[_tfreq] = dict({'fft':FT,'fft_freq':FQ-24,
                                'TOD':dat_i,'final_value':thisDat['%.1f'%_tfreq],
                                'measured_value':thisMeas['%.1f'%_tfreq],
                                'value_error':thisErr['%.1f'%_tfreq],
                                'valid':valid
                                })
        
    return thisDat,thisErr,outDict


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def full_model(w,tau,gamma,tau_0,Rnow,inductance):
    tau_e = 2.0*inductance/(Rnow)
#    model = lambda w,gamma,tau,tau_0,tau_e: abs((1./np.sqrt(1+w**2*tau_e**2))*(1./(((1+imag*w*tau_0)*(1+imag*w*tau_e))/(1+imag*w*tau_0/(1+gamma))+(tau_0/tau -1)*(1-imag*2*tau_e))))
    model = lambda w,tau,gamma,tau_0,tau_e: abs((1./np.sqrt(1+w**2*tau_e**2))*(1./(((1+imag*w*tau_0)*(1+imag*w*tau_e))/(1+(imag*w*tau_0)/(1+gamma))+(tau_0/tau-1)*(1-imag*2*tau_e))))
    norm_model = model(w,tau,gamma,tau_0,tau_e) / model(w[0],tau,gamma,tau_0,tau_e)
    return norm_model

def full_model_mp(params,w,data,err,fjac=None):
    tau,gamma,tau_0,Rnow,inductance = params
#    if tau_0 >= tau:
    output = (data-full_model(w,tau,gamma,tau_0,Rnow,inductance))/err
#    else:
#        output=np.ones(len(w))*1E10
    return 0, output
    
def fit_data(fit_type,freqs,amps,errs,
             f_3db=10.0,
             gamma = 300,
             tau_0 = 60E-3,
             Rnow = 0,
             inductance = 0,
             tau_0_range = None):
    ''' fit the specified model

    inputs: 
        freqs: a list of frequencies 
        amps: a list of amplitudes corresponding to the freqs
        tau_0: a guess at the tau of the pole you are fitting
    outputs:
        the fit tau
    '''
    w = 2*pl.pi*pl.asarray(freqs) # convert frequency into angular frequency 
    tau = 1/(float(f_3db)*pl.pi*2)

    functkw = {'w':w,'data':amps,'err':errs}
    # setup the parinfo so we can freeze parameters of the fit
    #        the limited=[1,0] with limits=[0,0] below makes all parameters positive
    parinfo = [{'value':0., 'fixed':0, 'limited':[1,0], 'limits':[0.,0.]} for i in range(5)] 
    parinfo[0]['parname'] = 'tau'
    parinfo[1]['parname'] = 'gamma'
    parinfo[2]['parname'] = 'tau_0'
    parinfo[3]['parname'] = 'Rnow'
    parinfo[4]['parname'] = 'inductance'

    # restrict gamma > 1
    parinfo[1]['limited'] = [1,0]
    parinfo[1]['limits'] = [1,0]

    if fit_type == 'single_body':
        parinfo[1]['fixed'] = 1 # freezes gamma
        parinfo[3]['fixed'] = 1 # freezes Rnow
        parinfo[4]['fixed'] = 1 # freezes inductance
        gamma = 1E10 # takes gamma->infinity
        Rnow = 1E10 # takes Rn->infinity
        inductance = 1E-10 # takes inductance->0
        parinfo[2]['tied'] = 'p[0]' # makes tau=tau_0

    elif fit_type == 'single_body_lr':
        parinfo[1]['fixed'] = 1 # freezes gamma
        parinfo[3]['fixed'] = 1 # freezes Rnow
        parinfo[4]['fixed'] = 1 # freezes inductance
        gamma = 1E10 # takes gamma->infinity
        parinfo[2]['tied'] = 'p[0]' # makes tau=tau_0
    elif fit_type == 'two_body':
        parinfo[3]['fixed'] = 1 # freezes Rnow
        parinfo[4]['fixed'] = 1 # freezes inductance
        Rnow = 1E10 # takes Rn->infinity
        inductance = 1E-10 # takes inductance->0

    elif fit_type == 'two_body_lr':
        parinfo[3]['fixed'] = 1 # freezes Rnow
        parinfo[4]['fixed'] = 1 # freezes inductance

        if tau_0_range:
            parinfo[2]['limited'] = [1,1]
            parinfo[2]['limits'] = tau_0_range
            tau_0 = np.average(tau_0_range)
        # restrict gamma 
        parinfo[1]['limited'] = [1,0]
        parinfo[1]['limits'] = [1,0]


    else:
        print 'You did not supply a correct model option to the fitter'
    
    p0 = [tau,gamma,tau_0,Rnow,inductance]        
    for i in range(len(p0)): parinfo[i]['value']=p0[i] 
#    print asdf
    fit_output = mpfit.mpfit(full_model_mp, p0, parinfo=parinfo,
                             functkw = functkw,quiet=True)
    tau,gamma,tau_0,Rnow,inductance =  fit_output.params
    chi2_reduced = fit_output.fnorm/fit_output.dof
    print fit_output.print_results()
    if fit_output.status<0:
        success = False
        tau_err,gamma_err,tau_0_err,Rnow_err,inductance_err =  [0,0,0,0,0]
    else:
        success = True
        tau_err,gamma_err,tau_0_err,Rnow_err,inductance_err =  fit_output.perror#*np.sqrt(chi_reduced)
        
    return success,tau, tau_err, gamma, gamma_err, tau_0, tau_0_err,Rnow,Rnow_err,inductance,inductance_err,chi2_reduced

def make_plots(data,savePlots=True,directory='',
               note='',
               fit_model = 'single_body',
               fit_freq_high_cut = None,
               fit_freq_low_cut = None,
               inductance = 0,
               Rt = 0,
               gamma = 200,
               tau_0 = 60E-3,
               tau_0_range = None,
               plot_residual = False,
               verbose=False):
    '''
    In addition to making visually pleasing plots, 
       this function calls the fitting routine for the data.

       plot_residual is not implemented yet!!!

    '''

    for i,bolo in enumerate(data.keys()):
        _now = data[bolo]
        #initiate a figure to use
        pl.figure(i+1,figsize=(10,10));pl.clf();
        if plot_residual:
            pl.subplot(121)
        pl.plot(1,1,'k',linestyle='none',label= '   %s   %s   %s         %s              %s       %s '%(r'$\frac{R_{TES}}{Rn}$',r'$f_{3db}\,(Hz)$',r'$\tau\,(\pm)\,[ms]$',r'$\gamma\,(\pm)$',r'$\tau_{o},(\pm)\,[ms]$',r'$\frac{\chi^{2}}{DOF}$'))
        freq_min_to_plot = 100
        freq_max_to_plot = 100
        min_amp_to_plot = 0.9

        pl.axhline(pl.sqrt(0.5),color='k',linestyle='--')#,alpha=0.7)

        # cycle through the responses at different Rfractions
        rfrac_count = len(_now.keys())
        for _nrfrac, rfrac in enumerate(sorted(_now.keys(),reverse=True)):
            freqs,responses,errors,valids = [],[],[],[]
            for _freq in np.sort(_now[rfrac]['tickle_frequencies'].keys()):
                _now_tick = _now[rfrac]['tickle_frequencies'][_freq]
                freqs.append(float(_freq))
                responses.append(_now_tick['final_value'])
                errors.append(_now_tick['value_error'])
                valids.append(_now_tick['valid'])
            # reorder the frequencies and responses to increasing in freq
#            freqs,responses,errors = pl.array(zip(*sorted(zip(freqs,responses,errors)))) # I dont think I need this if I order the keys above
                
            # only fit up to a specified frequency if fit_freq_high_cut and/or fit_freq_low_cut is supplied
            if fit_freq_high_cut:
                print 'Only fitting frequencies below %s.'%(str(fit_freq_high_cut))
                ind_to_fit = pl.find(np.array(freqs) < fit_freq_high_cut)[-1]
                pl.axvline(fit_freq_high_cut,color='k',alpha=0.3)#linestyle='--'
            else:
                ind_to_fit = -1
            if fit_freq_low_cut:
                print 'Only fitting frequencies above %s.'%(str(fit_freq_low_cut))
                start_ind = pl.find(np.array(freqs) >= fit_freq_low_cut)[0]
            else:
                start_ind = 0

            # normalize the response to the first VALID point
            ind_valid = valids[start_ind:].index(True)# this is the first 'valid' point
            start_ind += ind_valid
#            start_ind = valids.index(True)# this is the first 'valid' point
            respnorm = pl.asarray(responses)/responses[start_ind]
            errnorm = pl.asarray(errors)/responses[start_ind]
            freqs = pl.asarray(freqs)

            # finding the 3db point? T-funk, check this out more
            ii = min([pl.find(respnorm>0.5)[-1], len(respnorm)-2])    ### check this
            f_3db = pl.interp(pl.sqrt(0.5),respnorm[[ii+1,ii]],freqs[[ii+1,ii]])  ### check this

            # fit the data
            print '---------------------------------Now fitting rfrac = '+str(rfrac)
            fit_output = fit_data(fit_model,freqs[start_ind:ind_to_fit],
                                  respnorm[start_ind:ind_to_fit],
                                  errnorm[start_ind:ind_to_fit], f_3db = f_3db,
                                  tau_0 = tau_0,
                                  Rnow = Rt*float(rfrac), 
                                  inductance = inductance,
                                  tau_0_range = tau_0_range)

            fit_success,tau,tau_err,gamma,gamma_err,tau_0,tau_0_err,Rnow,Rnow_err,inductance,inductance_err,chi2_reduced = fit_output
            
            print 'tau = %6.3f +/- %6.8f ms'%(tau*1000,tau_err*1000)
            print 'tau_0 = %6.3f +/- %6.8f ms'%(tau_0*1000,tau_0_err*1000)
            print 'gamma = %6.2f +/- %6.2f'%(gamma,gamma_err)
            print 'Rnow = %6.4f +/- %6.5f Ohms'%(Rnow,Rnow_err)
            print 'Inductance = %6.2f +/- %6.5f uH'%(inductance*1E6,inductance_err*1E6)
            print 'chi^2_reduced = %6.4f'%(chi2_reduced)

            # if tau is unresponably high (corresponds to a timeconstant >500)
            #    set tau = NA in the label
            if 1/(2*pl.pi*tau) < 500 and fit_success==True:
                label = '%8.3f  %6.1f  %6.1f (%.2f)  %6.1f (%.2f)   %6.1f (%.2f)  %6.2f'%(float(rfrac),1/(2*pl.pi*tau),tau*1000.,tau_err*1000.,gamma,gamma_err,tau_0*1000.,tau_0_err*1000.,chi2_reduced)
            else:
                label = '%8.3f  %s  %6.1f   %s'%(float(rfrac),'  NA  ',tau*1000.,    ' fit invalid ')
#            print label,'+/-', tau_err*1000

            # plot the data
            pl.errorbar(freqs,respnorm,yerr=errnorm,fmt = '.',label=label,
                        color = pl.cm.jet(float(_nrfrac)/rfrac_count))
#                        color = pl.cm.gist_rainbow(float(_nrfrac)/rfrac_count))

            # make some nicely spaced points for the model to plot over
            freqs_plot = pl.logspace(1,np.log(freqs[-1])/np.log(freqs[start_ind]),base = freqs[start_ind])
            # plot the model 
            model_to_plot = full_model(2*pl.pi*freqs_plot,tau,gamma,tau_0,Rnow,inductance)
            if fit_success:
                pl.loglog(freqs_plot,model_to_plot,color=pl.gca().lines[-1].get_color(),linewidth=0.4)#,label = label)
                if plot_residual:
                    pl.subplot(122)
                    pl.loglog(freqs_plot,(model_to_plot-respnorm)/model_to_plot,color=pl.gca().lines[-1].get_color(),linewidth=0.4)#,label = label)
                    pl.subplot(121)
            if freq_min_to_plot>freqs[0]: freq_min_to_plot = freqs[0]
            if freq_max_to_plot<freqs[-1]: freq_max_to_plot = freqs[-1]
            if min_amp_to_plot>min(respnorm[np.isfinite(respnorm)]): min_amp_to_plot = min(respnorm[np.isfinite(respnorm)])
            
        pl.yscale('log'); pl.xscale('log')
        pl.legend(loc=3)
        pl.xlabel('freq (Hz)')
        pl.ylabel('Response (normalized to lowest freq)')
        pl.suptitle('%s'%directory)
        pl.title('Quick Look ETF (%s),  L=%.2guH  Rn=%.fmOhms'%(bolo,inductance*1E6,Rt*1E3))
        pl.xlim(freq_min_to_plot,freq_max_to_plot);pl.ylim(min_amp_to_plot,2)
        pl.grid()
        if savePlots:
            save_name = '%s/TauETF_MPIT_err_%s_%s.png'%(directory,bolo,note)
            pl.savefig(save_name)
            saved_pngs.append(save_name)
#            pl.savefig('%s/TauETF_MPIT_err_%s_%s'%(directory,bolo,note))
#            saved_pngs.append('%s/TauETF_MPIT_err_%s_%note.png'%(directory,bolo,note))



  
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

def analyze_etf_dir(directory, save=False, note='',plotFfts=False,
                    fit_model = 'single_body',
                    fit_freq_high_cut = None,
                    fit_freq_low_cut = None,
#                    Rn = None,
                    inductance = 20E-6,
                    gamma = 200,
                    tau_0 = 70E-3,
                    Rt = 0.261,
                    tau_0_range = None):
    '''  analyzes all of the etf data in a given directory

    fir_model options are:
         ['single_body', 'single_body_lr', 'two_body', 'two_body_lr']
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

        iv_pkls = glob('%s/%s_TuneBoloBombDan*.pkl'%(directory,_SQ))
        iv_times = [sptTime.toSptDatetime(_log.rpartition('Dan_')[-1].rpartition('.log')[0]).fil 
                      for _log in iv_pkls]
            
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
                data_points,err,outDict = get_etf_data_points(thisFL[ch],_SQ,ch,plot=plotFfts)
                if plotFfts:
                    raw_input('press any key (this pause is because plotFfts=True')
                try:
                    complete_crap = totalDict[_SQ+ch]
                except KeyError:
                    totalDict[_SQ+ch] = dict({rFrac[ch]:dict()})
                totalDict[_SQ+ch][rFrac[ch]]={'tickle_frequencies':outDict}
                data.setdefault(_SQ+ch,{})
                data[_SQ+ch][rFrac[ch]] = data_points

        make_plots(totalDict,savePlots=save,directory=directory,
                   note = note, 
                   fit_model=fit_model,
                   fit_freq_high_cut=fit_freq_high_cut,
                   fit_freq_low_cut=fit_freq_low_cut,
                   Rt=Rt,
                   inductance = inductance,
                   gamma = gamma,
                   tau_0_range = tau_0_range)
        if save:
            pkl.dump(totalDict,open( directory + '/TauETF_all_'+note+'.pkl','w'))
    
    print 'Sucessfully saved %i pngs' %(len(saved_pngs))
    for _png in saved_pngs:
        print _png
                                










##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
        



def investigate(data, savePlots=False,SqCh=None,Rfrac=None,
                fit_freq_high_cut = 1000,extras = [],checkTau = False,
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
            #   right now we are hard coding it to only fit with freqs below fit_freq_high_cut
                print 'Only fitting frequencies below %s'%(str(fit_freq_high_cut))
                ind_to_fit = pl.find(freqs < fit_freq_high_cut)[-1]
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

































































































########################################################################################
########################################################################################
#### Old models are below, they are correct, but do not use the full two body lr model
########################################################################################
########################################################################################

# def single_body_model(w,tau):#,norm_freq=None):
#     '''
#     This model will be normalized to the first given frequency
#     '''
#     return abs((1+imag*w[0]*tau)/(1+imag*w*tau))

# def single_body_model_mp(tau,w,data,err,fjac=None):
#     output = (data - single_body_model(w,tau))/err
#     status = 0
#     return status, output

# def body1_lr_model(w, tau_o, Rnow, inductance):
#     tau_e = 2.0*inductance/(Rnow)
#     model = 1./(np.sqrt(1+w**2*tau_e**2)) * (2*imag*tau_e*tau/tau_o + imag*w*tau - w**2*tau_e*tau + 1 - imag*w*tau_e)**(-1)
#     return model


# def two_body_model(w,tau,gamma,tau_0):
#     model = abs(1.0/( (1+imag*w*tau_0)/(1+imag*w*tau_0/(1+gamma)) + tau_0/tau -1))
#     norm_model = abs(model*( (1+imag*w[0]*tau_0)/(1+imag*w[0]*tau_0/(1+gamma)) + tau_0/tau -1))
#     return norm_model

# def two_body_model_mp(params,w,data,err,fjac=None):
#     tau,gamma,tau_0 = params
#     output = (data - two_body_model(w,tau,gamma,tau_0))/err
#     return 0, output
