from crython.cryoboard import CryoBoard, CryoStreamer
import pylab as pl
import numpy as np
import math,scipy.stats
from scipy.interpolate import UnivariateSpline
from glob import glob
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from pywtl.core.wtl_ConvertUtils import convert_squid
import sptpol_software.util.files as files

## Upgrade from rt_plot.  Reads in DfMUX and cryoboard parsers, fits alpha, loopgain, L(P_opt), etc. 

parser_dir = '20150606_cooldown/LC017/RT1'

def convert_temp(downsample=50,save=True):
  cb = CryoBoard('192.168.1.129')

  raw = glob('/home/cryo/Data/output/'+parser_dir+'/temp_'+str(downsample/25)+'sec'+'.txt')
  
  if len(raw) == 0:
    temp_v=pl.fromfile('/home/cryo/Data/output/'+parser_dir+'/cryoboard/voltage_rms_10.dat',np.float64)
  
    temp_i=pl.fromfile('/home/cryo/Data/output/'+parser_dir+'/cryoboard/current_rms_10.dat',np.float64)

    pad_size=math.ceil(float(temp_v.size)/downsample)*downsample-temp_v.size  ## Padding with zeros
    temp_v=np.append(temp_v,np.zeros(pad_size)*np.NaN)
    temp_i=np.append(temp_i,np.zeros(pad_size)*np.NaN)

    v=scipy.stats.nanmean(temp_v.reshape(-1,downsample), axis=1)  ## Downsampling
    i=scipy.stats.nanmean(temp_i.reshape(-1,downsample), axis=1)

    v=v*1.1957222e-9 # ADC to volts conversion
    i=i*1.1849500e-13 # ADC to amps conversion
    temp=v/i

    for _t,_r in enumerate(temp):
     temp[_t] = cb.convert_sensor_adc_resistance(11,_r,cb.OHMS,cb.KELVIN)
    
    if save==True:
      file=open('/home/cryo/Data/output/'+parser_dir+'/temp_'+str(downsample/25)+'sec'+'.txt','w')
      for item in temp:
        file.write("%s\n" % item)

  else:
    temp=pl.array(files.readColumns(raw[0],header=['temp'],column_header=True)['temp'])
  
  return temp

def convert_current(sq,ch='01',downsample=380): 

  if sq=='Sq1':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m1_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m1_c'+ch+'_q'
  if sq=='Sq2':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m2_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m2_c'+ch+'_q'
  if sq=='Sq3':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m3_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m3_c'+ch+'_q'
  if sq=='Sq4':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m4_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/139/192_168_1_139/192_168_1_139_bolo_m4_c'+ch+'_q'
  if sq=='Sq5':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m1_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m1_c'+ch+'_q'
  if sq=='Sq6':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m2_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m2_c'+ch+'_q'
  if sq=='Sq7':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m3_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m3_c'+ch+'_q'
  if sq=='Sq8':
    signal_i = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m4_c'+ch+'_i'
    signal_q = '/home/cryo/Data/output/'+parser_dir+'/137/192_168_1_137/192_168_1_137_bolo_m4_c'+ch+'_q'

  ii = pl.fromfile(signal_i,np.int32)
  iq = pl.fromfile(signal_q,np.int32)

  pad_size=math.ceil(float(ii.size)/downsample)*downsample-ii.size  ## Padding with zeros
  ii=np.append(ii,np.zeros(pad_size)*np.NaN)
  iq=np.append(iq,np.zeros(pad_size)*np.NaN)

  ii=scipy.stats.nanmean(ii.reshape(-1,downsample), axis=1)  ## Downsampling
  iq=scipy.stats.nanmean(iq.reshape(-1,downsample), axis=1)

  magnitude=ii
  for _i in range(len(ii)):
     magnitude[_i] = pl.sqrt(float(ii[_i])**2+float(iq[_i])**2)

  current = magnitude*convert_squid(units='Normalized').theo_AnToIbolo(1,2,N=16,eff=1)/2.0**(24-1)

  return current

def rt_plot(sqlst,chlst=['01'],avg_time=100,plot_alpha=False,s=0.01,plot_loopgain=False,fig=2,sweep=None,save=True):
  #Changed this to make plotting all R(T) curves easier when we're working with many detectors. 2015-02-17 DPD

  if sweep=='cooling':
    marker='v'
  elif sweep=='warming':
    marker='^'
  else:
    marker='v'

  if not isinstance(sqlst,list):
    sqlst = [sqlst]
  if not isinstance(chlst,list):
    chlst = [chlst]

  temp = convert_temp(downsample=avg_time*25,save=save)

  for num,sq in enumerate(sqlst):
    i=0.0
    pl.figure(fig+num)			# Plotting curves by squid
    pl.xlabel('Temperature (K)')
    pl.ylabel('R (Ohms)')
    pl.grid(b=True)
    pl.title(sq+' R-T curves')
    for ch in chlst:   
      current = convert_current(sq,ch=ch,downsample=avg_time*190)
      length_diff=len(current)-len(temp)
      resistance = ((.2641e-6)/current[length_diff:]) #.2641e-6
      pl.scatter(temp,resistance,c=pl.cm.jet(float(i/len(chlst))),label=sq+'-'+ch+'-W',marker=marker,edgecolors='none',s=40)
      i=i+1
    pl.legend(loc=4)


######### Fitting alpha #########
  if plot_alpha:
	    
    x=pl.sort(temp)*1e3
    y=pl.sort(resistance)

    ##### Use only the data within the transition #####
    ind=pl.where((x > 410) & (x < 440))
    x=x[ind]
    y=y[ind] 

    pl.figure(8)
    pl.clf()
    pl.scatter(pl.log(x),pl.log(y),s=2)

    fit = UnivariateSpline(pl.log(x),pl.log(y),s=s)
    xfit = pl.linspace (pl.log(x)[0],pl.log(x)[-1],1000)
    yfit=fit(xfit)
    
    pl.plot(xfit,yfit,'--')

    pl.figure(fig)
    pl.subplot(2,1,1)
    pl.scatter(x,pl.array(y)/pl.sort(y)[-1],s=3)#,color=colors[int(sq)],label='Sq'+sq)
    pl.plot(pl.e**pl.array(xfit),pl.e**pl.array(yfit)/pl.e**pl.array(yfit)[-1],'b-',label='Cooling')
    pl.xlim(410,440)
    pl.ylabel('R/Rn')
    pl.suptitle(sq+'_Ch'+ch+'RT measurement')
    pl.legend(loc=4)
    pl.grid(b=True)

    alpha=[]
    xfittemp=[]
    for _x in xfit:
      alpha.append(fit.derivatives(_x)[1])
      xfittemp.append(pl.e**_x)
    alpha=pl.array(alpha)
    xfittemp=pl.array(xfittemp)

    pl.subplot(2,1,2)
    pl.plot(xfittemp,alpha,'-b')
    pl.xlabel('Temperature (mK)')
    pl.ylabel('Alpha')
    pl.xlim(410,440)
    pl.grid(b=True)

######## Loopgain ###########
  if plot_loopgain:
    pl.figure(fig)
    pl.clf()

    if sq=='Sq2': # k and n come from G(T) measurements.
      k=322.      # P is transition bias power with stage at 315 mK
      n=2.7
      Rn=1.350
      P=44.5
    if sq=='Sq4':
      k=318.
      n=2.7
      Rn=1.350
      P=43.9
    if sq=='Sq5':
      k=267.
      n=2.7
      Rn=1.325
      P=37.5
    if sq=='Sq7':
      k=284.
      n=2.6
      Rn=1.296
      P=41.2
    if sq=='Sq8':
      k=284.
      n=2.6
      Rn=1.227
      P=39.0

   
    # iv_dir = '/home/cryo/Data/output/20140728_204918_IV_full/'
    # item=dict()

    # pkls = glob(iv_dir+sq+'_TuneBoloCombDan*.pkl')
    # if len(pkls)==0:
    #         print "No matching .pkl file found."
    #         return	
    # for _pkl in pkls:
    #         _label= _pkl.partition('IV_')[-1].rpartition('/')[0]
    #         item[_label] = [pkl.load(open(_pkl,'r'))]	

    # iv_rfracs=[]
    # P=[]	
    # for _run in item:
    #         _it = item[_run][0]
    #         vals = iv.analyzeIVData(_it['data'][1],iv_type='TuneBoloCombDan', 
    #                                 parRes = 0.08, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
    #         # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
    #         #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

    #         #ind=pl.argmin(vals['i'])
    #         iv_rfracs.append(vals['r_frac'])
    #         P.append(vals['p'])
    # iv_rfracs=pl.sort(pl.array(iv_rfracs)[0])
    # P=pl.sort(pl.array(P)[0])

    # Pbias=pl.interp(pl.sort(pl.e**pl.array(yfit)/pl.e**pl.array(yfit)[-1]),iv_rfracs,P)

    # #P = (0.23811**2)/(Rn*pl.e**pl.array(yfit)/pl.e**pl.array(yfit)[-1]) # Power in pW
    # loopgain=(alpha*Pbias)/(n*k*(xfittemp*1e-3)**n)

    # #pl.plot(pl.e**pl.array(yfit)/pl.e**pl.array(yfit)[-1],loopgain,'-r',label='Cooling')
    # pl.plot(Pbias+5.5,loopgain,'-b',label='Cooling')
  
    # pl.grid(b=True)
    # pl.xlabel('Pbias [pW]')
    # pl.ylabel('Loopgain')
    # #pl.title(sq+' 0.06 pW bias at 1 ohm')
    # #pl.title(sq+' Predicted Loopgain with Pbias='+str(P)+' pW')
    # pl.title(sq+' Loopgain with stage at 315mK')
    # pl.legend()
  
    ## Plotting Loopgain as function of optical power ##

    Psat = 9
    r_fracs=pl.array([0.5,0.6,0.7,0.8,0.9])
    Popt=pl.linspace(0,8.99,100)

    for i in range(len(r_fracs)):
      vbias2 = (Psat-6)*r_fracs[i]*Rn
      R = vbias2/(Psat-Popt)
      rf = (R/Rn)[pl.where(R/Rn < 1.0)[0]]
      a=pl.interp(rf,pl.e**pl.array(yfit)/pl.e**pl.array(yfit)[-1],alpha)
      loopgain=a*(Psat-Popt[0:len(rf)])/(n*Psat)
      pl.plot(Popt[0:len(rf)],loopgain,label='Rfrac='+str(r_fracs[i]))

    pl.xlabel('Total Optical Power [pW]')
    pl.ylabel('Loopgain')
    pl.title(sq+', biased with Popt = 6 pW')
    pl.grid(b=True)
    pl.legend()
