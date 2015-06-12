from crython.cryoboard import CryoBoard, CryoStreamer
import pylab as pl
import numpy as np
import math,scipy.stats
from scipy.interpolate import UnivariateSpline
from glob import glob
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from pywtl.core.wtl_ConvertUtils import convert_squid

## Upgrade from rt_plot.  Reads in DfMUX and cryoboard parsers, fits alpha, loopgain, L(P_opt), etc. 

parser_dir = '20141009_optical_cooldown/Hglamp_05cms'

def convert_temp(downsample=50):
  cb = CryoBoard('192.168.1.129')

  temp_v=pl.fromfile('/home/cryo/Data/output/'+parser_dir+'/cryoboard/voltage_rms_12.dat',np.float64)
  
  temp_i=pl.fromfile('/home/cryo/Data/output/'+parser_dir+'/cryoboard/current_rms_12.dat',np.float64)

  pad_size=math.ceil(float(temp_v.size)/downsample)*downsample-temp_v.size  ## Padding with zeros
  temp_v=np.append(temp_v,np.zeros(pad_size)*np.NaN)
  temp_i=np.append(temp_i,np.zeros(pad_size)*np.NaN)

  v=scipy.stats.nanmean(temp_v.reshape(-1,downsample), axis=1)  ## Downsampling
  i=scipy.stats.nanmean(temp_i.reshape(-1,downsample), axis=1)

  v=v*1.1957222e-9 # ADC to volts conversion
  i=i*1.1849500e-13 # ADC to amps conversion
  temp=v/i

  for _t,_r in enumerate(temp):
   temp[_t] = cb.convert_sensor_adc_resistance(13,_r,cb.OHMS,cb.KELVIN)

  return temp

def convert_current(sq):

  if sq=='Sq1':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m1_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m1_c01_q'
  if sq=='Sq2':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m2_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m2_c01_q'
  if sq=='Sq3':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m3_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m3_c01_q'
  if sq=='Sq4':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m4_c01_q'
  if sq=='Sq5':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m1_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m1_c01_q'
  if sq=='Sq6':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m2_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m2_c01_q'
  if sq=='Sq7':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m3_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m3_c01_q'
  if sq=='Sq8':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m4_c01_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m4_c01_q'

  ii = pl.fromfile(signal_i,np.int32)
  iq = pl.fromfile(signal_q,np.int32)


  magnitude=ii
  for _i in range(len(ii)):
     magnitude[_i] = pl.sqrt(float(ii[_i])**2+float(iq[_i])**2)

  current = magnitude*convert_squid(units='Normalized').theo_AnToIbolo(1,2,N=16,eff=1)/2.0**(24-1)

  return current*1e6

current = convert_current('Sq8')
pl.plot(current)
pl.ylabel('Current, muA')
pl.show()

