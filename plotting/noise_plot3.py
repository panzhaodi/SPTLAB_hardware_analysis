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
from crython.cryoboard import CryoBoard, CryoStreamer
import math,scipy.stats
from scipy.interpolate import UnivariateSpline
from glob import glob
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from pywtl.core.wtl_ConvertUtils import convert_squid

#parser_dir = 'SqCon027_SqBdB12/noise5'
#parser_dir = '20150327_cooldown/noise'
parser_dir = 'SqCon009_SqBdB15/noise5'
def noise_plot(sq,ch='01', r_frac='0.7', biasv=1.9):
  if sq=='Sq1':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m1_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m1_c'+ch+'_q'
  if sq=='Sq2':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m2_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m2_c'+ch+'_q'
  if sq=='Sq3':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m3_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m3_c'+ch+'_q'
  if sq=='Sq4':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m4_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139'+'/192_168_1_139/192_168_1_139_bolo_m4_c'+ch+'_q'
  if sq=='Sq5':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m1_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m1_c'+ch+'_q'
  if sq=='Sq6':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m2_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m2_c'+ch+'_q'
  if sq=='Sq7':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m3_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m3_c'+ch+'_q'
  if sq=='Sq8':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m4_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137'+'/192_168_1_137/192_168_1_137_bolo_m4_c'+ch+'_q'
  biasv=biasv*8.967529/1.9 #convert from bias amplitude to muV

  data_i= fromfile(signal_i,np.int32)[0:200000]
  data_q= fromfile(signal_q,np.int32)[0:200000]
  magnitude = fromfile(signal_i,np.int32)[0:200000]
  for i in range(0,len(data_i)):
     magnitude[i] = sqrt(abs(float(data_i[i])*float(data_i[i]))+abs(float(data_q[i])*float(data_q[i])))
  magnitude = (1e6)*magnitude*convert_squid(units='Normalized').theo_AnToIbolo(1,2,N=16,eff=1)/2.0**(24-1)
  sample_time=0.00524288
  n=len(data_i)
  time=np.linspace(0, (n-1)*sample_time, n)
  pl.figure(1)
  plot(time, magnitude)
  pl.xlabel('Time(s)')
  pl.ylabel('Current(uA)')
  pl.title(sq+' time ordered data'+' with V_bias=%.3f muV'%biasv)#+', r_frac='+r_frac)
  results={}
#  frequencies= np.fft.fftfreq(n, sample_time)[1:n/2]
  T=sample_time * n
  sampling =1.0/sample_time
  df =1.0/T
  spectrum= np.fft.fft(magnitude)
  if np.mod(n,2)==0:
     freq=np.array(range(0,n/2+1))/T
     psd_temp= np.zeros(len(freq))
     psd_temp[0]= abs(spectrum[0]*np.conj(spectrum[0]))
     psd_temp[n/2]= abs(spectrum[n/2]*np.conj(spectrum[n/2]))
     psd_temp[1:n/2]=abs(2*spectrum[1:n/2]*np.conj(spectrum[1:n/2]))
  else :
     freq=np.array(range(0,(n-1)/2+1))/T
     psd_temp= np.zeros(len(freq))
     psd_temp[0]= abs(spectrum[0]*np.conj(spectrum[0]))
     psd_temp[1:n/2+1]=abs(2*spectrum[1:n/2+1]*np.conj(spectrum[1:n/2+1]))
  norm=n*n*df
  psd_norm=psd_temp/norm

  frequencies=freq[1:n/2]
  spectrum2 = np.sqrt(8./3.)*1e6*np.sqrt(psd_norm)[1:n/2] #hanning window normalization factor, convert to pA
                       
  pl.figure()
  plt.loglog(frequencies, spectrum2)
  plt.title('Noise Spectrum for '+sq+' with V_bias=%.3f muV'%biasv)
  plt.xlabel('Frequency (Hz)')
  plt.ylabel('pA/rtHz'  )

  def fit_function(x, a,  floor ):
	    return np.sqrt(a*a/x/x+floor*floor)

  parameters2=curve_fit(fit_function,frequencies[10:len(frequencies)],spectrum2[10:len(frequencies)], p0=[4, 20] ) #Neglect the first 2 points, which will give us a bad fit.
  [a , floor]=parameters2[0]			
  fit=fit_function(frequencies,a, floor)


  pl.plot(frequencies, fit)
  pl.figtext(.20,.25,'Model = $\sqrt{(Af^{-1})^2+Floor^2}$',fontsize=14)
  pl.figtext(.20,.20,'A = %10.5f' % a,fontsize=14)
  pl.figtext(.20,.15,'Floor = %10.2f' % floor,fontsize=14)
  print np.average(spectrum2)


  '''
  def fit_function(x, a, n, floor ):
	    return np.sqrt(a*(x**n)*a*(x**n)+floor*floor)

  parameters2=curve_fit(fit_function,frequencies[10:len(frequencies)],spectrum2[10:len(frequencies)], p0=[4,-0.7, 20] ) #Neglect the first 2 points, which will give us a bad fit.
  [a, n , floor]=parameters2[0]			
  fit=fit_function(frequencies,a,n, floor)
  pl.plot(frequencies, fit)
  pl.figtext(.20,.30,'Model = $\sqrt{(Af^n)^2+Floor^2}$',fontsize=14)
  pl.figtext(.20,.25,'A = %10.5f' % a,fontsize=14)
  pl.figtext(.20,.20,'n = %10.2f' % n,fontsize=14)
  pl.figtext(.20,.15,'Floor = %10.2f' % floor,fontsize=14)
  print np.average(spectrum2)
  '''