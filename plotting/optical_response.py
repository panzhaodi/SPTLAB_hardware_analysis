import pylab as pl
import numpy as np
import cPickle as pkl
from scipy.interpolate import UnivariateSpline
from pylab import fromfile
from math import sqrt
from pylab import plot
import matplotlib.pyplot as plt
from pylab import loglog
from math import pow
from scipy.optimize import curve_fit
from glob import glob
from scipy.integrate import quad
#from lmfit import minimize, Parameters, Parameter, report_fit

def plot(rfrac,directory = '/home/cryo/Data/output/20141203_cooldown'):
	'''This program make the optical response of the detector as a function of chopper frequency. It is very similar to ETF measurement. An example call would be: optical_response.plot(0.8)
          '''
	item=dict()

	parseris = glob(directory+'/optical_response'+'/'+rfrac+'/*'+'/192_168_1_139/192_168_1_139_bolo_m2_c01_i')

	if len(parseris)==0:
	    print "No matching parser i file found."
	    return	
	for _parseri in parseris:
#look up the name of the folder
#The name of the folder is like /20140822_cooldown/optical_response/0.8/15/192_168_1_137/192_168_bolo_m4_c01_i, the label will be 0.8/15/i. 0.8 is the depth in transition, 15 is the chopper frequency, i is the type of signal we are looking at.
	    _label= _parseri.partition('optical_response/')[-1].partition('/192_168_1_139')[0]+'/i'
            item[_label] = fromfile(_parseri,np.int32)	
#?does this need to be a list


        parserqs = glob(directory+'/optical_response'+'/'+rfrac+'/*'+'/192_168_1_139/192_168_1_139_bolo_m2_c01_q')

	if len(parserqs)==0:
	    print "No matching parser q file found."
	    return	
	for _parserq in parserqs:
#look up the name of the folder
#The name of the folder is like /20140822_cooldown/optical_response/0.8/15/192_168_bolo_m4_c01_i, the label will be 0.8/15/i. 0.8 is the depth in transition, 15 is the chopper frequency, i is the type of signal we are looking at.
	    _label= _parserq.partition('optical_response/')[-1].partition('/192_168_1_139')[0]+'/q'
            item[_label] = fromfile(_parserq,np.int32)	
#?does this need to be a list
#find the python script to find the maximum member's index
# You will need to change the frequency list if you have made any change to it.
#	freq_list=['5','10','20','30','40','50','60','70','80','90','100','120','150']
        freq_list=['10','20','30','40','50','60','70','85','90']
        conversion_factor=1.932651e-5
        accurate_freq=[]
        amplitude_list=[]
        freq_error=[]
        amp_error=[]
	for fr in freq_list:	
	   data_i= item[rfrac+'/'+fr+'/i']
           data_q= item[rfrac+'/'+fr+'/q']
	   magnitude1 = np.linspace(0,0,15440)
           for i in range(0,15440):
               magnitude1[i] = sqrt(abs(float(data_i[i])*float(data_i[i]))+abs(float(data_q[i])*float(data_q[i])))
           magnitude1 = conversion_factor*magnitude1
           sample_time=0.00524288
           sample_number=15440
           time=np.linspace(0, (sample_number-1)*sample_time, sample_number)
           frequencies= np.fft.fftfreq(sample_number, sample_time)[5:sample_number/2]
           spectrum1 = abs(np.fft.fft(magnitude1*np.hanning(sample_number)))[5:sample_number/2] # ignore DC
           spectrum1 = spectrum1 * (1e6 # uA -> pA
                               * np.sqrt(2.) # definition of ASD, see below
                               * np.sqrt(8./3.) # renormalize to compensate for hanning window
                               / np.sqrt(sample_number) # amplitude spectrum convention
                               / np.sqrt(190.734863281/2)) # go to amplitude spectral density                         

           spectrum1= list(spectrum1)
           ind1= spectrum1.index(max(spectrum1))
           amp1 = np.sum(spectrum1[ind1-2: ind1+3])
#find the index of the max element. denoted by ind1
           frequencies=list(frequencies)
           frequency1=frequencies[ind1]

           if fr=='90':
            spectrum1= list(spectrum1)[7000:]
            ind1= spectrum1.index(max(spectrum1))
            amp1 = np.sum(spectrum1[ind1-2: ind1+3])
#find the index of the max element. denoted by ind1
            frequencies=list(frequencies)[7000:]
            frequency1=frequencies[ind1]             


#calculate for the second segment

	   magnitude2 = np.linspace(0,0,15440)
           for i in range(0,15440):
               magnitude2[i] = sqrt(abs(float(data_i[i+15440])*float(data_i[i+15440]))+abs(float(data_q[i+15440])*float(data_q[i+15440])))
           magnitude2 = conversion_factor*magnitude2
           sample_time=0.00524288
           sample_number=15440
           time=np.linspace(0, (sample_number-1)*sample_time, sample_number)
           frequencies= np.fft.fftfreq(sample_number, sample_time)[5:sample_number/2]
           spectrum2 = abs(np.fft.fft(magnitude2*np.hanning(sample_number)))[5:sample_number/2] # ignore DC
           spectrum2 = spectrum2 * (1e6 # uA -> pA
                               * np.sqrt(2.) # definition of ASD, see below
                               * np.sqrt(8./3.) # renormalize to compensate for hanning window
                               / np.sqrt(sample_number) # amplitude spectrum convention
                               / np.sqrt(190.734863281/2)) # go to amplitude spectral density                         
           spectrum2= list(spectrum2)
           ind2= spectrum2.index(max(spectrum2))
#find the index of the max element. denoted by ind1
           amp2 = np.sum(spectrum2[ind2-2: ind2+3])
           frequencies=list(frequencies)
           frequency2=frequencies[ind2]
           if fr=='90':
            spectrum2= list(spectrum2)[7000:]
            ind2= spectrum2.index(max(spectrum2))
            amp2 = np.sum(spectrum2[ind2-2: ind2+3])
#find the index of the max element. denoted by ind1
            frequencies=list(frequencies)[7000:]
            frequency2=frequencies[ind2]  
#calculate for the third segment
	   magnitude3 = np.linspace(0,0,15440)
           for i in range(0,15440):
               magnitude3[i] = sqrt(abs(float(data_i[i+30880])*float(data_i[i+30880]))+abs(float(data_q[i+30880])*float(data_q[i+30880])))
           magnitude3 = conversion_factor*magnitude3
           sample_time=0.00524288
           sample_number=15440
           time=np.linspace(0, (sample_number-1)*sample_time, sample_number)
           frequencies= np.fft.fftfreq(sample_number, sample_time)[5:sample_number/2]
           spectrum3 = abs(np.fft.fft(magnitude3*np.hanning(sample_number)))[5:sample_number/2] # ignore DC
           spectrum3 = spectrum3 * (1e6 # uA -> pA
                               * np.sqrt(2.) # definition of ASD, see below
                               * np.sqrt(8./3.) # renormalize to compensate for hanning window
                               / np.sqrt(sample_number) # amplitude spectrum convention
                               / np.sqrt(190.734863281/2)) # go to amplitude spectral density                         
#find the index of the max element. denoted by ind3
           spectrum3= list(spectrum3)
           ind3= spectrum3.index(max(spectrum3))
           amp3 = np.sum(spectrum3[ind3-2: ind3+3])
           frequencies=list(frequencies)
           frequency3=frequencies[ind3]

           if fr=='90':
            spectrum3= list(spectrum3)[7000:]
            ind3= spectrum3.index(max(spectrum3))
            amp3 = np.sum(spectrum3[ind3-2: ind3+3])
#find the index of the max element. denoted by ind1
            frequencies=list(frequencies)[7000:]
            frequency3=frequencies[ind3] 

           accurate_freq.append(np.average([frequency1, frequency2, frequency3]))
           amplitude_list.append(np.average([amp1,amp2,amp3]))
           freq_err= np.std([frequency1, frequency2, frequency3])/sqrt(2.0)
           amp_err= np.std([amp1,amp2,amp3])/sqrt(2.0)
           freq_error.append(freq_err)
           amp_error.append(amp_err)
           pl.plot(frequencies, spectrum1)
        rescale=amplitude_list[0]
        amplitude_list=np.array(amplitude_list)/rescale
        amp_error=np.array(amp_error)/rescale
       # pl.errorbar(accurate_freq, amplitude_list, xerr=freq_error, yerr=amp_error, fmt='o')
        pl.xlabel('Chopper frequency(Hz)')
        pl.ylabel('Amplitude(arbitrary units)')
 #axis
        pl.title('Optical response at '+rfrac+'R_normal')

'''
        def fit_function(x, t0,A ):
	    return A/np.sqrt(1.0+4*3.14159*3.1415926*t0**2*np.array(x)**2)
        parameters=curve_fit(fit_function,accurate_freq,amplitude_list, p0=[0.01,1000] ) 
        [t0,A]=parameters[0]
        frequencies=np.linspace(3,92,1000)
        fit=fit_function(frequencies,t0,A)

  #      pl.plot(frequencies, fit)


        params = Parameters()
        params.add('t0', value= 0.01, min=0.0001, max=1)
        params.add('A', value= 30000 , min=10000, max=40000)

        x= np.array(accurate_freq)
        amp = np.array(amplitude_list)
        def fit_function(params, x, amp):
            """ model decaying sine wave, subtract data"""
            t0 = params['t0'].value
            A = params['A'].value
            return  amp-A/np.sqrt(1.0+4*3.14159*3.1415926*t0**2*np.array(x)**2)



# do fit, here with leastsq model
        result = minimize(fit_function, params, args=(x,amp))
# calculate final result
        def fit_function(x, t0,A ):
	    return A/np.sqrt(1.0+4*3.14159*3.1415926*t0**2*np.array(x)**2)
        [t0,A]=[params['t0'].value, params['A'].value]
        frequencies=np.linspace(3,92,1000)
        fit=fit_function(frequencies,t0,A)
# write error report
        pl.plot(frequencies, fit, color=colors, label=rfrac+' R_normal  '+'A = %2.0f, $t_0$=%2.4f s'% (A,t0))
        pl.legend()
        report_fit(params)
# try to plot results
#pl.plot(x, final, 'r',label='Fit curve')


        pl.figtext(.20,.20,'Model = $A/\sqrt{1+\omega^2 t_0^2}$',fontsize=20)
  
 # pl.figtext(.20,.20,'n = %10.2f' % n,fontsize=14)
 # pl.figtext(.20,.15,'Floor = %10.2f' % floor,fontsize=14)

     '''      





