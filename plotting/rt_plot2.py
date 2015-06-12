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

#default directory
directory = '20150606_cooldown/RT5'
#output directory
output_dir='20150606_cooldown'

def convert_temp(downsample=50,save=True,parser_dir=directory):
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

def convert_current(sq,ch='01',downsample=190.734863281,parser_dir=directory): 

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
  ii_d=[]
  iq_d=[]
  for i in range(int(len(ii)/downsample)):
     start=int(i*downsample)
     end=int((i+1)*downsample)
   
     ii_d.append(np.mean(ii[start:end]))
     iq_d.append(np.mean(iq[start:end]))
  ii=np.array(ii_d)
  iq=np.array(iq_d)
  magnitude=ii
  for _i in range(len(ii)):
     magnitude[_i] = pl.sqrt(float(ii[_i])**2+float(iq[_i])**2)

  current = magnitude*convert_squid(units='Normalized').theo_AnToIbolo(1,2,N=16,eff=1)/2.0**(24-1)

  return current

def rt_plot(sqlst,chlst=['01'],avg_time=100,plot_alpha=False,s=0.01,plot_loopgain=False,fig=2,save=True, warming_folder=None, cooling_folder=None):
  #Changed this to make plotting all R(T) curves easier when we're working with many detectors. 2015-02-17 DPD
  '''
  This function is for plotting the RT curve of a detector. 
  Args:
        sqlst : the list of SQUIDs we are trying to plot. Example: ['Sq1', 'Sq2']
	chlist : the list of channel in each SQUID we want to plot. Example ['01']
	avg_time : the downsample time
	plot_alpha : the option to plot alpha or not
	s : the interpolation parameter
	plot_loopgain : the option to plot loopgain or not
	sweep : whether you are taking the data during warming or cooling
	save : whether you want to save the temperature to save time for the future
	warming_folder : the directory where you store the warming parser data
	cooling_folder : the directory where you store the cooling parser data
  '''
  if not isinstance(sqlst,list):
    sqlst = [sqlst]
  if not isinstance(chlst,list):
    chlst = [chlst]
  color_list= [ 'chocolate', 'hotpink','gold', 'red', 'blue', 'green','orchid', 'burlywood','darksalmon', 'gray','darkviolet', 'yellowgreen', 'whitesmoke']

  if warming_folder!=None and cooling_folder!=None:
       for num,sq in enumerate(sqlst):
    	  i=0.0
    	  pl.figure(fig+num)			# Plotting curves by squid
    	  pl.xlabel('Temperature (K)')
    	  pl.ylabel('R (Ohms)')
    	  pl.grid(b=True)
    	  pl.title(sq+' R-T curves')
          temp_w = convert_temp(downsample=avg_time*25,save=save, parser_dir=output_dir+'/'+warming_folder)
          temp_c = convert_temp(downsample=avg_time*25,save=save, parser_dir=output_dir+'/'+cooling_folder)
    	  for ch in chlst:   
      		current_w = convert_current(sq,ch=ch,downsample=avg_time*190.734863281, parser_dir=output_dir+'/'+warming_folder)
                resistance_w=((.2641e-6)/current_w)
                if len(resistance_w)>len(temp_w):
                      length_diff=len(current_w)-len(temp_w)
      		      resistance_w = resistance_w[:len(resistance_w)-length_diff] 
                else: 
                      length_diff=len(temp_w)-len(current_w)
                      temp_w = temp_w[:len(temp_w)-length_diff] 
                current_c = convert_current(sq,ch=ch,downsample=avg_time*190.734863281, parser_dir=output_dir+'/'+cooling_folder)
                resistance_c=((.2641e-6)/current_c)
                if len(resistance_c)>len(temp_c):
                      length_diff=len(current_c)-len(temp_c)
      		      resistance_c = resistance_c[:len(resistance_c)-length_diff] 
                else: 
                      length_diff=len(temp_c)-len(current_c)
                      temp_c = temp_c[:len(temp_c)-length_diff] 


      		pl.scatter(temp_w,resistance_w,c=color_list[int(i)],label=sq+'Ch'+ch,marker='o',edgecolors='none',s=40)
                pl.scatter(temp_c,resistance_c,c='k',marker='.',s=40)
                pl.figtext(0.6,0.25,'Colored dots - warming', fontsize=18)
                pl.figtext(0.6,0.22, 'Black dots - cooling', fontsize=18)
      		i=i+1
    	  pl.legend(loc=4)

  elif ((warming_folder!=None and cooling_folder==None) or (warming_folder==None and cooling_folder!=None)):
       if warming_folder!=None and cooling_folder ==None:
          folder_label='(warming)'
          parser_d=output_dir+'/'+warming_folder
       else:
          folder_label='(cooling)' 
          parser_d=output_dir+'/'+cooling_folder

       for num,sq in enumerate(sqlst):
    	  i=0.0
    	  pl.figure(fig+num)			# Plotting curves by squid
    	  pl.xlabel('Temperature (K)')
    	  pl.ylabel('R (Ohms)')
    	  pl.grid(b=True)
    	  pl.title(sq+' R-T curves'+folder_label)
          temp = convert_temp(downsample=avg_time*25,save=save, parser_dir=parser_d)

    	  for ch in chlst:   
      		current = convert_current(sq,ch=ch,downsample=avg_time*190.734863281, parser_dir=parser_d)
                resistance=((.2641e-6)/current)
                if len(resistance)>len(temp):
                      length_diff=len(current)-len(temp)
      		      resistance = resistance[length_diff:] 
                else: 
                      length_diff=len(temp)-len(current)
                      resistance = resistance[length_diff:] 
              

      		pl.scatter(temp,resistance,c=color_list[int(i)],label=sq+'Ch'+ch,marker='o',edgecolors='none',s=40)


      		i=i+1
    	  pl.legend(loc=4)
  else:
       print "No parser directory given."


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
