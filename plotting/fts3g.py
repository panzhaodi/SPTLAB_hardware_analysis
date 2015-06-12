import pylab as pl
import sptpol_software.util.files as files
import cPickle as pkl
from glob import glob
from pyfts_package import *
import numpy as np

# For chopped scans, you might have to shift the frequencies by 60 Hz. I'm not sure why this is yet.

def fts(parser_dir,sq,speed=1.0,chop=False,save=False,return_data=False):
   
  data_dir='/home/cryo/Data/output/20141203_cooldown/'
  rate = 191.   #sampling rate in Hz
  #speed = optical velocity in cm/s

  band = [70,110]   #approximate band of the detector, used to find interferrorgrams in code

  if sq=='Sq1':
    signal_i = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m1_c01_i'
    signal_q = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m1_c01_q'
  if sq=='Sq2':
    signal_i = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m2_c01_i'
    signal_q = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m2_c01_q'
  if sq=='Sq3':
    signal_i = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m3_c01_i'
    signal_q = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m3_c01_q'
  if sq=='Sq4':
    signal_i = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m4_c01_i'
    signal_q = data_dir+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m4_c01_q'
  if sq=='Sq5':
    signal_i = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m1_c01_i'
    signal_q = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m1_c01_q'
  if sq=='Sq6':
    signal_i = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m2_c01_i'
    signal_q = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m2_c01_q'
  if sq=='Sq7':
    signal_i = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m3_c01_i'
    signal_q = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m3_c01_q'
  if sq=='Sq8':
    signal_i = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m4_c01_i'
    signal_q = data_dir+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m4_c01_q'
  
  raw = glob(data_dir+parser_dir+'/*'+sq+'.txt')
  if len(raw) == 0:
    data_i= pl.fromfile(signal_i,np.int32)
    data_q= pl.fromfile(signal_q,np.int32)
    
    magnitude=data_i
    for _i in range(len(data_i)):
      magnitude[_i] = pl.sqrt(float(data_i[_i])**2+float(data_q[_i])**2)

    if save:
      file=open(data_dir+parser_dir+'/'+parser_dir.rpartition('/')[-1]+'_'+sq+'.txt','w')
      for item in magnitude:
         file.write("%s\n" % item)
  
  else:
    magnitude=files.readColumns(raw[0],header=['magnitude'],column_header=True)['magnitude']
    
  cts=pl.array(magnitude)
  if chop:
    cts=-1*cts
  #pl.figure()
  #pl.clf()
  #pl.plot(range(len(cts)),cts)   #plot it
  #cts = cts[1.5e4:*]  # if you want to cut some data points off the beginning


  out = analyze_fts_scan(cts, speed, band=band, rate=rate, chop=chop,hann=True, pband = [10, 400],absv=False)  #process the FTS data# if it's high enough signal-to-noise drop the /  abs at the end

  avg = add_fts_scan(out)  #average the fts scans

  pl.figure()
  pl.clf()
  pl.plot(avg['freq'], avg['real']/max(avg['real'])) # plot the average
  pl.xlim(10, 400)
  pl.xlabel('Frequency (GHz)')
  pl.ylabel('Normalized Response')
  pl.plot(avg['freq'], avg['im']/max(avg['real']), 'b--') #overplot the imagninary component to get an estimate of the noise
  pl.grid(b=True)
  pl.title('Averaged Spectra')

  if return_data:
    if save:
      pkl.dump(avg,open(data_dir+parser_dir+'/'+parser_dir.rpartition('/')[-1]+'_'+sq+'_spectrum.pkl','wb'))
    return avg
  #look at write_png.pro write_jpg.pro (write_jpeg?.pro)
  #the colors with my color table can look really weird. 

  #openw, 1, 'Hglamp_140GHz_1cms_spectrum.txt'
  #for i=0, 1910 do printf, 1,string(avg.freq[i]), string(avg.real[i]), string(avg.im[i])
  #close,1

