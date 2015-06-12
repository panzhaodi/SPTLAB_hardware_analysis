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
#sample call:  noise_plot.noise_plot('2','0.7')
def noise_plot(sq='2', r_frac='0.7'):
  if sq =='2' and r_frac =='0.7':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13977/192_168_1_139/192_168_1_139_bolo_m2_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13977/192_168_1_139/192_168_1_139_bolo_m2_c01_q'
  if sq=='4' and r_frac=='0.7':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13977/192_168_1_139/192_168_1_139_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13977/192_168_1_139/192_168_1_139_bolo_m4_c01_q'
  if sq=='5' and r_frac=='0.7':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13777/192_168_1_137/192_168_1_137_bolo_m1_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13777/192_168_1_137/192_168_1_137_bolo_m1_c01_q'
  if sq=='7' and r_frac=='0.7':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13777/192_168_1_137/192_168_1_137_bolo_m3_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13777/192_168_1_137/192_168_1_137_bolo_m3_c01_q'
  if sq=='8' and r_frac=='0.7':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13777/192_168_1_137/192_168_1_137_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13777/192_168_1_137/192_168_1_137_bolo_m4_c01_q'



  if sq =='2' and r_frac =='0.6':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13966/192_168_1_139/192_168_1_139_bolo_m2_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13966/192_168_1_139/192_168_1_139_bolo_m2_c01_q'
  if sq=='4' and r_frac=='0.6':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13966/192_168_1_139/192_168_1_139_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13966/192_168_1_139/192_168_1_139_bolo_m4_c01_q'
  if sq=='5' and r_frac=='0.6':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13766/192_168_1_137/192_168_1_137_bolo_m1_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13766/192_168_1_137/192_168_1_137_bolo_m1_c01_q'
  if sq=='7' and r_frac=='0.6':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13766/192_168_1_137/192_168_1_137_bolo_m3_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13766/192_168_1_137/192_168_1_137_bolo_m3_c01_q'
  if sq=='8' and r_frac=='0.6':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13766/192_168_1_137/192_168_1_137_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13766/192_168_1_137/192_168_1_137_bolo_m4_c01_q'





  if sq =='2' and r_frac =='0.5':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13955/192_168_1_139/192_168_1_139_bolo_m2_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13955/192_168_1_139/192_168_1_139_bolo_m2_c01_q'
  if sq=='4' and r_frac=='0.5':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13955/192_168_1_139/192_168_1_139_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13955/192_168_1_139/192_168_1_139_bolo_m4_c01_q'
  if sq=='5' and r_frac=='0.5':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13755/192_168_1_137/192_168_1_137_bolo_m1_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13755/192_168_1_137/192_168_1_137_bolo_m1_c01_q'
  if sq=='7' and r_frac=='0.5':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13755/192_168_1_137/192_168_1_137_bolo_m3_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13755/192_168_1_137/192_168_1_137_bolo_m3_c01_q'
  if sq=='8' and r_frac=='0.5':
    signal_i = '/home/cryo/Data/output/20140723_cooldown/13755/192_168_1_137/192_168_1_137_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/20140723_cooldown/13755/192_168_1_137/192_168_1_137_bolo_m4_c01_q'
 
  conversion_factor=1.932651e-5
  data_i= fromfile(signal_i,np.int32)
  data_q= fromfile(signal_q,np.int32)
  magnitude = fromfile(signal_i,np.int32)
  for i in range(0,len(data_i)):
     magnitude[i] = sqrt(abs(float(data_i[i])*float(data_i[i]))+abs(float(data_q[i])*float(data_q[i])))
  magnitude = conversion_factor*magnitude
  sample_time=0.00524288
  sample_number=len(data_i)
  time=np.linspace(0, (sample_number-1)*sample_time, sample_number)
  pl.figure(1)
  plot(time, magnitude)
  pl.xlabel('Time(s)')
  pl.ylabel('Current(uA)')
  pl.title(' Sq'+sq+' time ordered data'+', r_frac='+r_frac)
  results={}
  frequencies= np.fft.fftfreq(sample_number, sample_time)[3:sample_number/2]
  spectrum = abs(np.fft.fft(magnitude*np.hanning(sample_number)))[3:sample_number/2] # ignore DC
  spectrum = spectrum * (1e6 # uA -> pA
                         * np.sqrt(2.) # definition of ASD, see below
                         * np.sqrt(8./3.) # renormalize to compensate for hanning window
                         / np.sqrt(sample_number) # amplitude spectrum convention
                         / np.sqrt(190.734863281/2)) # go to amplitude spectral density
#  results['spectrum'] = spectrum                           
  pl.figure(2)
  plt.loglog(frequencies, spectrum)
  plt.title('Noise Spectrum for Sq'+sq+' at r_frac='+r_frac)
  plt.xlabel('Frequency (Hz)')
  plt.ylabel('pA/rtHz'  )

  def fit_function(x, a, n, floor ):
	    return np.sqrt(a*(x**n)*a*(x**n)+floor*floor)

  parameters2=curve_fit(fit_function,frequencies[3:len(frequencies)],spectrum[3:len(frequencies)], p0=[4,-0.7, 20] ) #Neglect the first 2 points, which will give us a bad fit.
  [a, n , floor]=parameters2[0]			
  fit=fit_function(frequencies,a,n, floor)
  loglog(frequencies, fit)
  pl.figtext(.20,.30,'Model = $\sqrt{(Af^n)^2+Floor^2}$',fontsize=14)
  pl.figtext(.20,.25,'A = %10.5f' % a,fontsize=14)
  pl.figtext(.20,.20,'n = %10.2f' % n,fontsize=14)
  pl.figtext(.20,.15,'Floor = %10.2f' % floor,fontsize=14)



