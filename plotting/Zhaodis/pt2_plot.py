'''
Module that contains two different functions for plotting Power vs Temperature data.
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
import pylab as pl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from math import pow
from scipy.integrate import simps

def plot_stage_pt(sq,chan=1,temp='???',parRes=0.1,legends=False):
	'''Example call: pt_plot.plot_stage_pt('Sq4',2,'300')
	This function plots Pbias at turn-around at various bath temperatures and fits it to P = k*(Tc**n-Tb**n)'''

	path = '/home/cryo/Data/output/20140511_cooldown' # Directory with labeled data, e.g. /20140421_181415_IV_343mK
	item=dict()

	pkls = glob(path+'/*'+temp+'mK/'+sq+'_TuneBoloCombDan*.pkl')
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
    	pl.clf()
	pl.scatter(Tb,pTurn)
	pl.xlabel('Bath Temperature [K]')
        pl.ylabel('Power on TES [pW]')
	if legends:        
	    pl.legend()
	pl.title(sq +'_Ch'+str(chan))
	pl.grid()
	
	#def fit_function(x,a,b,c,d,e):
	    #return (a*x**4 + b*x**3 + c*x**2 + d*x + e)

	#def fit_function(x,a,b):
	    #return a*x+b

	def fit_function(Tb,k,Tc,n):
	    return k*(Tc**n-Tb**n)

	parameters=curve_fit(fit_function,Tb,pTurn,p0=[370,.515,3])
	[k,Tc,n]=parameters[0]	
	#[a,b,c,d,e]=parameters[0]
	#[a,b]=parameters[0]		
	x = pl.linspace(.250,.550,100)
	#fit = fit_function(x,a,b,c,d,e)
	#fit=fit_function(x,a,b)	
	fit=fit_function(x,k,Tc,n)
	
	pl.figtext(.15,.37,'$P=k({T_c}^n-{T_B}^n)$',fontsize=16)
	pl.figtext(.15,.30,'k = %.0f$\pm$ %.0f' % (k, pl.sqrt(parameters[1][0,0])),fontsize=14)
	pl.figtext(.15,.25,'Tc = %.3fK$\pm$ %.3fK' % (Tc, pl.sqrt(parameters[1][1,1])),fontsize=14)
	pl.figtext(.15,.20,'n = %.1f$\pm$ %.2f' % (n, pl.sqrt(parameters[1][2,2])),fontsize=14)
	pl.plot(x,fit)
	#print "[a,b,c,d,e]= ",[a,b,c,d,e]
	#print "[a,b]= ",[a,b]
	print "[k,Tc,n]= ",[k,Tc,n]	
	#print pl.sqrt(parameters[1][1,1])
	#print pl.sqrt(parameters[1][4,4])
	print pl.sqrt(parameters[1])
	#pl.show(block=False)

def plot_cl_pt(sq,chan=1,r_frac = .80, parRes = {'Sq3': .1,'Sq4':.1, 'Sq8':0.08},legends=True):
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots Pbias at a specified Rfrac with different coldload temps and fits it'''

	path = '/home/cryo/Data/output/20140723_cooldown/coldload' # Directory containing labeled data folders, e.g. .../20140426_162036_IV_08.51K
	item=dict()
	if not isinstance(sq,list):
	  sq=[sq]



 	with open("220_s21_filter.txt") as f:
 	   data=f.read()

 	data = data.split('\n')
 	x = [row.split(' ')[0] for row in data]
 	del data[-1]
 	x = [row.split(' ')[0] for row in data]
 	y=[row.split(' ')[-1].split('\r')[0] for row in data]
 	z = [row.split(' ')[1] for row in data]
 	for i in range (0,5001):
  	  x[i]= float(x[i])

 	for i in range (0,5001):
          y[i]= float(y[i])

 	for i in range (0,5001):
  	 z[i] = float (z[i])

 	w=[pow(float(a),2)+pow(float(b),2) for a,b in zip(y,z)]
 	filterf = interp1d(x, w)

 	with open("antenna_radiation_220.txt") as f:
    	 data=f.read()

 	data = data.split('\n')
 	x = [row.split(' ')[0] for row in data]
 	del data[-1]
 	x = [row.split(' ')[0] for row in data]
 	y=[row.split(' ')[1]  for row in data]
 
 	for i in range (0,1001):
  	 x[i]= float(x[i])

 	for i in range (0,1001):
  	 y[i]= 2*float(y[i])
 
 	antenna = interp1d(x, y)



 	with open("fts220sq6.txt") as f:
    	 data=f.read()

 	data = data.split('\n')
        del data[len(data)-1]
        data=[row.split('\r')[0] for row in data]
 	x = [row.split(',')[0] for row in data]
 	y = [row.split(',')[1] for row in data]
 	
 
 	for i in range (0,len(x)):
  	 x[i]= float(x[i])

 	for i in range (0,len(x)):
  	 y[i]= float(y[i])
 
 	fts = interp1d(x, y)




 	def fit_function(T,a,b):
   	 ans=[] 
   	 for _num in T:
      	   x=pl.linspace(140,300,10000)
      	   y= fts(x)*b*2*5e11*((6.626e-34)*x*1e9*1e9)/(pl.e**(((4.8e-11)*x*1e9)/_num)-1)
      	   ans.append(a-simps(y,x))	     
   	 return pl.array(ans)

	

#	    ans.append(a-b*5e11*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),200e9,240e9)[0])
#	    ans.append(a - (b**2)*5e11*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),200e9,240e9)[0]) - (c**2)*1e12*quad(lambda x:((6.626e-34)*(x**3)/((3e8)**2))/(pl.e**(((4.8e-11)*x)/_num)-1),350e9,500e9)[0])

#	colors=['b','g']
	colors=['b','g','r']
#	colors=['k','m']

	pl.figure(2)

        for i,j in enumerate(sq):
  	  pkls = glob(path+'/*??.??K/'+j+'_TuneBoloCombDan*.pkl')
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
                                    parRes = parRes[j], cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	      if vals['hasPTurn']:
		if float(_run.rpartition('_')[-1].partition('K')[0]) > 5: # In case you only want to plot certain temp ranges.
  	          p=pl.interp(r_frac,vals['r_frac'],vals['p'])
	          pR.append(p)
	          Tcl.append(float(_run.rpartition('_')[-1].partition('K')[0]))
		
	  pR=pl.array(pR)
	  Tcl=pl.array(Tcl)

 	  parameters=curve_fit(fit_function,Tcl,pR)

	  [a,b]=parameters[0]
#	  [a,b,c,n]=parameters[0]
#	  [a,c,n]=parameters[0]
			
	  x = pl.linspace(.1,61,1000)
	
	  fit=fit_function(x,a,b)
#	  fit =fit_function(x,a,b,c,n)
#	  fit=fit_function(x,a,c,n)
          pl.figure(1)
#	  pl.scatter(Tcl,pR,edgecolors=colors[i],facecolors='none',marker='o')#,label=j+': a = %.3f, b = %.3f'%(a,b))
#	  pl.scatter(Tcl,pR,color=colors[i],label=j+r': a = %.3f pW, b = %.3f, c = %.3f, n = %.3f'%(a,b**2,c**2,n))
	  pl.scatter(Tcl,pR,color=colors[i])#,label=j+': a = %.3f, c = %.3f, n = %.3f'%(a,c,n))


#	  red_chi_sq=(sum((pR-fit_function(Tcl,a,b))**2))/(len(pR)-2)		# These Chi-squareds don't make sense yet
#	  red_chi_sq=(sum((pR-fit_function(Tcl,a,b,c,n))**2))/(len(pR)-4)
#	  red_chi_sq=(sum((pR-fit_function(Tcl,a,c,n))**2))/(len(pR)-3)


	  pl.plot(x,fit,color=colors[i],linestyle='-',label = j+': a = %.4f,b = %.6f'%(a,b))
        #pl.figtext(.15,.17, r'$P_{bias} = a-b*10^{12}{0.95^4}\int^{\nu_1}_{\nu_0}FTS(\nu)\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu$',fontsize=18)
	pl.figtext(.15,.17, r'$P_{bias} = a-b*10^{12} {0.95^4} \int^{\nu_1}_{\nu_0}F(\nu)A(\nu)\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=22)
#	pl.figtext(.15,.30, r'$P_{bias} = a-10^{12}\frac{b}{2}\int^{\nu_1}_{\nu_0}\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu - cT^n$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=18)
#	pl.figtext(.15,.40, r'$P_{bias} = a - cT^n$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=18)
	pl.xlabel('Coldload Temperature [K]')
        pl.ylabel('Bias Power on TES [pW]')
	if legends:        
	    pl.legend(loc='best')
	pl.grid()
	x=pl.linspace(180,260,3000)
        pl.figure(2)
        pl.plot(x, 0.646663*fts(x))
        pl.xlabel('Frequency(GHz)')
        pl.ylabel('Efficiency')
        pl.title('Rescaled FTS data')
        pl.figure(3)
        pl.plot(x, antenna(x))
        pl.xlabel('Frequency(GHz)')
        pl.ylabel('Efficiency')
        pl.title('Antenna response(simulated by Gensheng)')
        pl.figure(4)
        pl.plot(x, filterf(x))
        pl.xlabel('Frequency(GHz)')
        pl.ylabel('Efficiency')
        pl.title('Filter function(simulated by Gensheng)')

        pl.figure(5)
        pl.plot(x, 0.646663*fts(x)/(0.95*0.95*0.95*0.95*filterf(x)*antenna(x)))
        pl.xlabel('Frequency(GHz)')
        pl.ylabel('Rescaled_FTS/(Antenna*Filter)')
        pl.title('Rescaled_FTS/(0.95^4*Antenna*Filter)')

        pl.figure(6)
        pl.plot(x,  (0.95*0.95*0.95*0.95*filterf(x)*antenna(x)))
        pl.xlabel('Frequency(GHz)')
        pl.ylabel('0.95^4*Antenna*Filter')
        pl.title('0.95^4*Antenna*Filter')

        x=pl.linspace(140,300,10000)
        opticalpower1= simps(0.646663*fts(x),x)
        opticalpower2= 0.97*0.95*0.95*0.95*0.95*simps(antenna(x)*filterf(x), x)
        print 'Normalized FTS integral is %10.6f, multiplication of antenna and filter is %10.6f ' %(opticalpower1, opticalpower2)




def plot_diff_powers(sq,parRes = {'Sq3': .1,'Sq4':.1},chan=1,legends=True):
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots the ratio of difference powers w.r.t 5 K at a specified Rfrac with different coldload temps'''

	path = '/home/cryo/Data/output/20140511_cooldown' # Directory containing labeled data folders, e.g. .../20140426_162036_IV_08.51K
	item=dict()
	data=dict()
	if not isinstance(sq,list):
	  sq=[sq]

	pl.figure(2)

        for i,j in enumerate(sq):
  	  pkls = glob(path+'/*??.??K/'+j+'_TuneBoloCombDan*.pkl')
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
                                    parRes = parRes[j], cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
            # values in vals: 'pturn','resp','loopgain','slopeSign','freq',
            #             'vRaw','i','r_frac','p','r','hasPTurn','v','iRaw'

            #ind=pl.argmin(vals['i'])
	      if vals['hasPTurn']:
		if float(_run.rpartition('_')[-1].partition('K')[0]) < 50: # In case you only want to plot certain temp ranges.
  	          ind=pl.where((vals['r_frac'] < .90) & (vals['r_frac'] > 0.75))
	          pR.append(vals['p'][ind].mean())
	          Tcl.append(float(_run.rpartition('_')[-1].partition('K')[0]))
		
	  pR=pl.array(pR)
	  print pR
	  Tcl=pl.array(Tcl)
	  T_min=pl.where(Tcl==Tcl.min())
	  print Tcl
	  data[j]=(pR-pR[T_min])
	  print data[j]

	
	pl.scatter(Tcl,(data['Sq4']/data['Sq3']),color='r')
	return data['Sq4']/data['Sq3']

def plot_stable_pt(sq,chan=1,r_frac = .80, parRes = {'Sq3': .1,'Sq4':.1},directory = '/home/cryo/Data/output/20140511_cooldown/ptime',legends=True,do_fit=True,cutoff=50):
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots Pbias at a specified Rfrac with different coldload temps and fits it. Assumes directory has labeled folders containing the IV data, e.g. .../20140426_162036_IV_08.51K'''

	
	item=dict()
	if not isinstance(sq,list):
	  sq=[sq]
		
 	with open("220_s21_filter.txt") as f:
 	   data=f.read()

 	data = data.split('\n')
 	x = [row.split(' ')[0] for row in data]
 	del data[-1]
 	x = [row.split(' ')[0] for row in data]
 	y=[row.split(' ')[-1].split('\r')[0] for row in data]
 	z = [row.split(' ')[1] for row in data]
 	for i in range (0,5001):
  	  x[i]= float(x[i])

 	for i in range (0,5001):
          y[i]= float(y[i])

 	for i in range (0,5001):
  	 z[i] = float (z[i])

 	w=[pow(float(a),2)+pow(float(b),2) for a,b in zip(y,z)]
 	filterf = interp1d(x, w)

 	with open("antenna_radiation_220.txt") as f:
    	 data=f.read()

 	data = data.split('\n')
 	x = [row.split(' ')[0] for row in data]
 	del data[-1]
 	x = [row.split(' ')[0] for row in data]
 	y=[row.split(' ')[1]  for row in data]
 
 	for i in range (0,1001):
  	 x[i]= float(x[i])

 	for i in range (0,1001):
  	 y[i]= 2*float(y[i])
 
 	antenna = interp1d(x, y)




 	def fit_function(T,a,b):
   	 ans=[] 
   	 for _num in T:
      	   x=pl.linspace(150,270,10000)
      	   y= b*0.95*0.95*0.95*0.95*antenna(x)*filterf(x)*5e11*((6.626e-34)*x*1e9*1e9)/(pl.e**(((4.8e-11)*x*1e9)/_num)-1)
      	   ans.append(a-simps(y,x))	     
   	 return pl.array(ans)
	colors=['b','g','r']

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
                                    parRes = parRes[j], cgain = _it['carrier_gain'], ngain = _it['nuller_gain'])
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

	  
	  facecolors = colors[i]#'none'
	  if do_fit:
    	      parameters=curve_fit(fit_function,Tcl,pR)
	      facecolors = colors[i]
     	      [a,b]=parameters[0]

	      x = pl.linspace(.1,40,1000)
	
	      fit=fit_function(x,a,b)

	      pl.plot(x,fit,color=colors[i],linestyle='-',label = j+': a = %.2f, b = %.3f $\pm$ %.3f'%(a,b,pl.sqrt(parameters[1][1,1])))

	  pl.scatter(Tcl,pR,color=colors[i],facecolor=facecolors)#,label=j+': a = %.3f, b = %.3f'%(a,b))


       
	pl.figtext(.15,.17, r'$P_{bias} = a-b*10^{12}\frac{0.95^4}{2}\int^{\nu_1}_{\nu_0}F(\nu)A(\nu)\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=22)
	pl.xlabel('Coldload Temperature [K]')
        pl.ylabel('Bias Power on TES [pW]')
	if legends:        
	    pl.legend()
	pl.grid(b=True)
