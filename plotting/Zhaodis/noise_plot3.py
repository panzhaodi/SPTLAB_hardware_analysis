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

parser_dir = '20150124_load_cooldown/noise_275mk'
#parser_dir = '20150124_load_cooldown/noise'

def noise_plot(sq='Sq1',ch='01', r_frac='0.7'):
  if sq=='Sq1':
    signal_i = '/home/cryo/Data/output/stage_test/200a200/192_168_1_139/192_168_1_139_bolo_m1_c01_i' 
    signal_q = '/home/cryo/Data/output/stage_test/200a200/192_168_1_139/192_168_1_139_bolo_m1_c01_q'


  data_i= fromfile(signal_i,np.int32)
  print len(data_i)
  data_q= fromfile(signal_q,np.int32)
  magnitude = fromfile(signal_i,np.int32)
  for i in range(0,len(data_i)):
     magnitude[i] = sqrt(abs(float(data_i[i])*float(data_i[i]))+abs(float(data_q[i])*float(data_q[i])))
  magnitude = (1e6)*magnitude*convert_squid(units='Normalized').theo_AnToIbolo(1,2,N=16,eff=1)/2.0**(24-1)
  sample_time=0.00524288
  n=len(data_i)

  time=np.linspace(0, (n-1)*sample_time, n)
  print n
  pl.figure(1)
  plot(time, magnitude)
  pl.xlabel('Time(s)')
  pl.ylabel('Current(uA)')
  pl.title(sq+' time ordered data'+', r_frac='+r_frac)

