'''
Module that contains various functions for plotting Power vs Temperature data.
plot_stage_pt plots Pbias at turn-around at various stage temperatures.
plot_cl_pt plots Pbias at a specified Rfrac with different coldload temperatures.
Both functions assume you have the data labeled in its own directory.
'''
import pylab as pl
import cPickle as pkl
import pywtl.common.analysis.an_iv as iv
from glob import glob
from scipy.optimize import curve_fit
from scipy.integrate import quad
import numpy as np
from scipy import interpolate

def plot_stage_pt(sq,chan=1,temp='???',directory = '/home/cryo/Data/output/20141230_cooldown',parRes=0.25,legends=False):
	'''Example call: pt_plot.plot_stage_pt('Sq4',2,'300')
	This function plots Pbias at turn-around at various bath temperatures and fits it to P = k*(Tc**n-Tb**n). Assumes directory has labeled folders containing the IV data, e.g. /20140421_181415_IV_343mK'''
	
	item=dict()

	pkls = glob(directory+'/*'+temp+'mK/'+sq+'_TuneBoloCombDan*.pkl')
	if len(pkls)==0:
	    print "No matching .pkl file found."
	    return	
	for _pkl in pkls:
	    _label= _pkl.partition('IV_')[-1].rpartition('/')[0]
            item[_label] = [pkl.load(open(_pkl,'r'))]	

	pTurn=[]
	Tb=[]	
	for _run in item:
            _it = item[_run][0]
            vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	    if vals['hasPTurn']:
	      pTurn.append(vals['pturn'])
	      Tb.append(int(_run.rpartition('_')[-1].partition('mK')[0])*10**(-3))
		
	pTurn=pl.array(pTurn)
	Tb=pl.array(Tb)
	pl.figure()
#    	pl.clf()
	pl.scatter(Tb,pTurn)
	pl.xlabel('Bath Temperature [K]')
        pl.ylabel('Power on TES [pW]')
	if legends:        
	    pl.legend()
	pl.title(sq +'_Ch'+str(chan))
	pl.grid(b=True)
	
	#def fit_function(x,a,b,c,d,e):
	    #return (a*x**4 + b*x**3 + c*x**2 + d*x + e)

	#def fit_function(x,a,b):
	    #return a*x+b

	def fit_function(Tb,k,Tc,n):
	    return k*(Tc**n-Tb**n)

	parameters=curve_fit(fit_function,Tb,pTurn,p0=[370,.430,3])
	[k,Tc,n]=parameters[0]	
	#[a,b,c,d,e]=parameters[0]
	#[a,b]=parameters[0]		
	x = pl.linspace(.250,.450,100)
	#fit = fit_function(x,a,b,c,d,e)
	#fit=fit_function(x,a,b)	
	fit=fit_function(x,k,Tc,n)
	
	pl.figtext(.15,.47,'$P=k({T_c}^n-{T_B}^n)$',fontsize=16)
	pl.figtext(.15,.40,'k = %.0f$\pm$ %.0f' % (k, pl.sqrt(parameters[1][0,0])),fontsize=14)
	pl.figtext(.15,.35,'Tc = %.3fK$\pm$ %.3fK' % (Tc, pl.sqrt(parameters[1][1,1])),fontsize=14)
	pl.figtext(.15,.30,'n = %.1f$\pm$ %.2f' % (n, pl.sqrt(parameters[1][2,2])),fontsize=14)
	pl.figtext(.15,.20,'G(Tc) = %.0f pW/K' % (n*k*(Tc**(n-1))),fontsize=14)
	pl.plot(x,fit)
	#print "[a,b,c,d,e]= ",[a,b,c,d,e]
	#print "[a,b]= ",[a,b]
	print "[k,Tc,n]= ",[k,Tc,n]	
	#print pl.sqrt(parameters[1][1,1])
	#print pl.sqrt(parameters[1][4,4])
	print pl.sqrt(parameters[1])
	return parameters[0]

def plot_cl_pt(sq,chan=1,r_frac = .80, parRes = .25,directory = '/home/cryo/Data/output/20140723_cooldown/all_coldload',legends=False,do_fit=False,cutoff=45):
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots Pbias at a specified Rfrac with different coldload temps and fits it. Assumes directory has labeled folders containing the IV data, e.g. .../20140426_162036_IV_08.51K'''

	
	item=dict()
	if not isinstance(sq,list):
	  sq=[sq]
		
	def fit_function(T,a,b):
          ans=[] 
	  for _num in T:
	    ans.append(a-b*1e12*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),200e9,240e9)[0])
#	    ans.append(a - (b**2)*5e11*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),200e9,240e9)[0]) - (c**2)*1e12*quad(lambda x:((6.626e-34)*(x**3)/((3e8)**2))/(pl.e**(((4.8e-11)*x)/_num)-1),350e9,500e9)[0])
#	    ans.append(a-c*(_num**n))
	  return pl.array(ans)

#	colors=['b','g']
	colors=['b','g','r', 'c', 'm', 'y', 'k', 'w']
#	colors=['k','m']

	pl.figure(2)

        for i,j in enumerate(sq):
  	  pkls = glob(directory+'/*??.??K/'+j+'_TuneBoloCombDan*.pkl')
  	  if len(pkls)==0:
  	      print "No matching .pkl file found."
  	      return
  	  for _pkl in pkls:
	      _label= _pkl.partition('IV_')[-1].rpartition('/')[0]
              item[_label] = [pkl.load(open(_pkl,'r'))]	

	  pR=[]
	  Tcl=[]	

	  for _run in item:
              _it = item[_run][0]
              vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	      if vals['hasPTurn']:
		if float(_run.rpartition('min')[-1].partition('K')[0]) < cutoff: # In case you only want to plot certain temp ranges.
  	          p=pl.interp(r_frac,vals['r_frac'],vals['p'])
	          pR.append(p)
	          Tcl.append(float(_run.rpartition('min')[-1].partition('K')[0]))
		
	  pR=pl.array(pR)
	  Tcl=pl.array(Tcl)

	  error=pl.sqrt((.01*pR)**2) # just some guess
	  
	  facecolors = 'none'
	  if do_fit:
    	      parameters=curve_fit(fit_function,Tcl,pR,sigma=error)
	      facecolors = colors[i]
     	      [a,b]=parameters[0]
#	      [a,b,c,n]=parameters[0]
#	      [a,c,n]=parameters[0]
			
	      x = pl.linspace(.1,45,1000)
	
	      fit=fit_function(x,a,b)
#	      fit =fit_function(x,a,b,c,n)
#	      fit=fit_function(x,a,c,n)

	      pl.plot(x,fit,color=colors[i],linestyle='-',label = j+': a = %.2f, b = %.2f $\pm$ %.2f'%(a,b,pl.sqrt(parameters[1][1,1])))

	  pl.scatter(Tcl,pR,color=colors[i],facecolor=facecolors)#,label=j+': a = %.3f, b = %.3f'%(a,b))
#	  pl.scatter(Tcl,pR,color=colors[i],label=j+r': a = %.3f pW, b = %.3f, c = %.3f, n = %.3f'%(a,b**2,c**2,n))
#	  pl.scatter(Tcl,pR,color=colors[i],label=j+': a = %.3f, c = %.3f, n = %.3f'%(a,c,n))


#	  red_chi_sq=(sum((pR-fit_function(Tcl,a,b))**2))/(len(pR)-2)		# These Chi-squareds don't make sense yet
#	  red_chi_sq=(sum((pR-fit_function(Tcl,a,b,c,n))**2))/(len(pR)-4)
#	  red_chi_sq=(sum((pR-fit_function(Tcl,a,c,n))**2))/(len(pR)-3)




	#pl.figtext(.15,.17, r'$P_{bias} = a-10^{12}\frac{b}{2}\int^{\nu_1}_{\nu_0}\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=22)
#	pl.figtext(.15,.30, r'$P_{bias} = a-10^{12}\frac{b}{2}\int^{\nu_1}_{\nu_0}\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu - cT^n$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=18)
#	pl.figtext(.15,.15, r'$P_{bias} = a - cT^n$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=18)
	pl.xlabel('Coldload Temperature [K]')
        pl.ylabel('Bias Power on TES [pW]')
	if legends:        
	    pl.legend()
	pl.grid(b=True)
	
def plot_diff_powers(sq,parRes = .25,chan=1,legends=True):           # For stripline loss measurements
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots the ratio of difference powers w.r.t 5 K at a specified Rfrac with different coldload temps'''

	path = '/home/cryo/Data/output/20140626_cooldown' # Directory containing labeled data folders, e.g. .../20140426_162036_IV_08.51K
	item=dict()
	data=dict()
	error=dict()
	if not isinstance(sq,list):
	  sq=[sq]

	pl.figure(2)

        for i,j in enumerate(sq):
  	  pkls = glob(path+'/*_??.??K/'+j+'_TuneBoloCombDan*.pkl')
  	  if len(pkls)==0:
  	      print "No matching .pkl file found."
  	      return
  	  for _pkl in pkls:
	      _label= _pkl.partition('IV_')[-1].rpartition('/')[0] # _label would be e.g. 05.45 K
              item[_label] = [pkl.load(open(_pkl,'r'))]	

	  pR=[]
	  Tcl=[]	

	  for _run in item:
              _it = item[_run][0]
              vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	      if vals['hasPTurn']:
		if float(_run.rpartition('min')[-1].partition('K')[0]) < 70: # In case you only want to plot certain temp ranges.
  	          ind=pl.where((vals['r_frac'] < .90) & (vals['r_frac'] > 0.75))
	          pR.append(vals['p'][ind].mean())
	          Tcl.append(float(_run.rpartition('_')[-1].partition('K')[0]))
		
	  pR=pl.array(pR)
	  print pR
	  Tcl=pl.array(Tcl)
	  T_min=pl.where(Tcl==Tcl.min())
	  print Tcl
	  data[j]=(pR-pR[T_min])
	  error[j] = pl.sqrt((.01*pR)**2+(.01*pR[T_min])**2)
	  print data[j]

	ind = data['Sq3'].nonzero()
	
	final_data = (data['Sq4'][ind])/(data['Sq3'][ind])
	final_error = pl.sqrt((((error['Sq4'][ind])/(data['Sq4'][ind]))**2 + ((error['Sq3'][ind])/(data['Sq3'][ind]))**2)*final_data**2)

	def fit_function(x,b):
	    return b

	parameters=curve_fit(fit_function,Tcl[ind],final_data,sigma=final_error)
	x=pl.linspace(0,60,2)	
	pl.plot(x,[parameters[0],parameters[0]],color='b',linewidth=2,label=' %.3f $\pm$ %.3f' %(parameters[0],pl.sqrt(parameters[1])))

	pl.errorbar(Tcl[ind],final_data,yerr=final_error,fmt='ko')
	pl.ylim(0,1.5)
	pl.xlabel('Coldload Temperature')
	pl.ylabel(r'$\frac{\Delta P_4}{\Delta P_3}$',fontsize=22,rotation='horizontal')
	pl.grid(b=True)
	pl.legend(loc='best')
	return Tcl[ind],final_data,parameters

def plot_cl_ptime(sq,chan=[1],temp=10,r_frac = .75, parRes = .25 ,directory = '/home/cryo/Data/output/20141230_cooldown',legends=True):
	'''Example call: pt_plot.plot_cl_ptime(['Sq1','Sq3'],chan=[1,2])
	This function plots Pbias at a specified Rfrac with different coldload temps. Assumes directory has labeled folders containing the IV data, e.g. .../20140426_162036_IV_60min40K'''

	if not isinstance(sq,list):
	  sq=[sq]

        if not isinstance(chan,list):
	  chan=[chan]

	colors=['b','g','r', 'c', 'm', 'y', 'k', 'w']

	pl.figure(1)

        for i,j in enumerate(sq):
	  item=dict()
  	  pkls = glob(directory+'/*min'+str(temp)+'*/'+j+'_TuneBoloCombDan*.pkl')
  	  if len(pkls)==0:
  	      print "No matching .pkl file found."
  	      return
  	  for _pkl in pkls:
	      _label= _pkl.partition('IV_')[-1].rpartition('/')[0]
              item[_label] = [pkl.load(open(_pkl,'r'))]	

	  pR=[]
	  time=[]	

	  for _run in item:
              _it = item[_run][0]
              vals = iv.analyzeIVData(_it['data'][chan[i]],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	      if vals['hasPTurn']:
  	          p=pl.interp(r_frac,vals['r_frac'],vals['p'])
	          pR.append(p)
	          time.append(float(_run.partition('min')[0]))
		  
	  pR=pl.array(pR)
	  time=pl.array(time)

	  pl.scatter(time,pR,color=colors[i],facecolor=colors[i],label=j)

	pl.xlabel('Minutes elapsed')
        pl.ylabel('Bias Power on TES [pW]')
	pl.title('Stage at 319mK, Coldload at %sK'%temp)
	if legends:        
	    pl.legend()
	pl.grid(b=True)

def plot_stable_pt(sq,chan=1,r_frac = .80, parRes = .25,directory = '/home/cryo/Data/output/20140626_cooldown/ptime',legends=True,do_fit=True,cutoff=50, use_sim=False, do_average=False,colorin='b'):
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots Pbias at a specified Rfrac with different coldload temps and fits it. Assumes directory has labeled folders containing the IV data, e.g. .../20140426_162036_IV_08.51K. Cutoff is the temperature up to which you want to plot. use_sim is to choose whether you want to use simulated filter function when fitting. do_average is to choose whether you want to average multiple measurements near a temperature'''

	
	
	if not isinstance(sq,list):
	  sq=[sq]
		
	#def fit_function(T,a,b):
        #  ans=[] 
	#  for _num in T:
	#    ans.append(a-b*1e12*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),100e9,350e9)[0])

	#  return pl.array(ans)
        #if we use Donna's simulation to fit the data. 
        # for 90GHz use donnasim90ghz.s2p, for 220GHz use donnasim220ghz.s2p
        
	colors=['b','g','r','k','m']

	pl.figure(2)

        for i,j in enumerate(sq):
	  item=dict()
  	  pkls = glob(directory+'/*??.??K/'+j+'_TuneBoloCombDan*.pkl')
  	  if len(pkls)==0:
  	      print "No matching .pkl file found."
  	      return
  	  for _pkl in pkls:
	      _label= _pkl.partition('IV_')[-1].rpartition('/')[0]
              item[_label] = [pkl.load(open(_pkl,'r'))]	

	  pR=[]
	  Tcl=[]	

	  for _run in item:
              _it = item[_run][0]
              vals = iv.analyzeIVData(_it['data'][chan],iv_type='TuneBoloCombDan', 
                                    parRes = parRes, cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	      if vals['hasPTurn']:
		if float(_run.rpartition('min')[-1].partition('K')[0]) < cutoff: # In case you only want to plot certain temp ranges.
  	          p=pl.interp(r_frac,vals['r_frac'],vals['p'])
	          pR.append(p)
	          Tcl.append(float(_run.rpartition('min')[-1].partition('K')[0]))
		
	  pR=pl.array(pR)
	  Tcl=pl.array(Tcl)
          #if you want to average multiple measurements near a temperature
	  if do_average:
             templist=[5.63,7.50,9.00,10.00,15.00,20.00,25.00,30.00,35.00,40.00]
             pRtemp=[]
             Tcltemp=[]
             Perror=[]
             Terror=[]
             for temp_item in templist:
                ptemp=[]
                ttemp=[]
                for u in range(0,len(Tcl)):
                #select the temperatures that are close to the temperature list for averaging
                  if abs(temp_item-Tcl[u])<0.1:
                    ptemp.append(pR[u])
                    ttemp.append(Tcl[u])
                #Put the averaged power and temperature in a list
                pRtemp.append(np.average(ptemp))
                Tcltemp.append(np.average(ttemp))
                #To calculate the error of the data. Now it may return error because some lists only have 1 element
                Perror.append(np.std(ptemp)/np.sqrt(len(ptemp)-1))
                Terror.append(np.std(ttemp)/np.sqrt(len(ttemp)-1))
             pR=pRtemp
             Tcl=Tcltemp
             

	  facecolors = 'none'
	  if do_fit:
             
    	      parameters=curve_fit(fit_function,Tcl,pR)#, sigma=Perror)
	      facecolors = colors[i]
     	      [a,b]=parameters[0]

	      x = pl.linspace(.1,40,1000)
	
	      fit=fit_function(x,a,b)

	      pl.plot(x,fit,color=colors[i],linestyle='-',label = j+': a = %.3f, b = %.3f $\pm$ %.3f'%(a,b,pl.sqrt(parameters[1][1,1])))

	  pl.scatter(Tcl,pR,color=colorin,facecolor=facecolors, label='Depth:%.3f'%r_frac, s=10, marker='.')#,label=j+': a = %.3f, b = %.3f'%(a,b))
          if do_average:
             pl.errorbar(Tcl, pR,  yerr=Perror, fmt='.')


	#pl.figtext(.15,.17, r'$P_{bias} = a-10^{12}b\int^{\nu_1}_{\nu_0}Filter$-$sim(\nu)\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=22)
	pl.xlabel('Coldload Temperature [K]')
        pl.ylabel('Bias Power on TES [pW]')
	if legends:        
	    pl.legend()
	pl.grid(b=True)