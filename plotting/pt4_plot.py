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


def plot_cl_pt(legends=True,do_fit=True,cutoff=42):
	'''Example call: pt_plot.plot_cl_pt(['Sq3','Sq4'])
	This function plots Pbias at a specified Rfrac with different coldload temps and fits it. Assumes directory has labeled folders containing the IV data, e.g. .../20140426_162036_IV_08.51K'''

	
		
	def fit_function(T,a,b):
          ans=[] 
	  for _num in T:
	    ans.append(a-b*5e11*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),137.5e9,160e9)[0])
#	    ans.append(a - (b**2)*5e11*quad(lambda x:((6.626e-34)*x)/(pl.e**(((4.8e-11)*x)/_num)-1),200e9,240e9)[0]) - (c**2)*1e12*quad(lambda x:((6.626e-34)*(x**3)/((3e8)**2))/(pl.e**(((4.8e-11)*x)/_num)-1),350e9,500e9)[0])
#	    ans.append(a-c*(_num**n))
	  return pl.array(ans)

#	colors=['b','g']
	colors=['r','g','b']
#	colors=['k','m']

	pl.figure(2)

	pR=[23.3485,23.1897,23.0160,22.7409,22.5227,22.2746, 22.0292]
	Tcl=[5.40,10.02,15.02,25.04,30.05,35.50, 40.00]	

 		
	pR=pl.array(pR)
	Tcl=pl.array(Tcl)

	error=pl.sqrt((.01*pR)**2) # just some guess
	  
	facecolors = 'none'
	if do_fit:
    	      parameters=curve_fit(fit_function,Tcl,pR,sigma=error)
	      facecolors = colors[0]
     	      [a,b]=parameters[0]
#	      [a,b,c,n]=parameters[0]
#	      [a,c,n]=parameters[0]
			
	      x = pl.linspace(.1,42,1000)
	
	      fit=fit_function(x,a,b)
#	      fit =fit_function(x,a,b,c,n)
#	      fit=fit_function(x,a,c,n)

	      pl.plot(x,fit,color=colors[0],linestyle='-',label = 'Sq3: a = %.2f, b = %.2f $\pm$ %.2f'%(a,b,pl.sqrt(parameters[1][1,1])))

	pl.scatter(Tcl,pR,color=colors[0],facecolor=facecolors)#,label=j+': a = %.3f, b = %.3f'%(a,b))
#	pl.scatter(Tcl,pR,color=colors[i],label=j+r': a = %.3f pW, b = %.3f, c = %.3f, n = %.3f'%(a,b**2,c**2,n))
#	pl.scatter(Tcl,pR,color=colors[i],label=j+': a = %.3f, c = %.3f, n = %.3f'%(a,c,n))


#	red_chi_sq=(sum((pR-fit_function(Tcl,a,b))**2))/(len(pR)-2)		# These Chi-squareds don't make sense yet
#	red_chi_sq=(sum((pR-fit_function(Tcl,a,b,c,n))**2))/(len(pR)-4)
#	red_chi_sq=(sum((pR-fit_function(Tcl,a,c,n))**2))/(len(pR)-3)




	pl.figtext(.15,.17, r'$P_{bias} = a-10^{12}\frac{b}{2}\int^{\nu_1}_{\nu_0}\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=22)
#	pl.figtext(.15,.30, r'$P_{bias} = a-10^{12}\frac{b}{2}\int^{\nu_1}_{\nu_0}\frac{h\nu}{e^{h\nu/{kT}}-1}d\nu - cT^n$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=18)
#	pl.figtext(.15,.15, r'$P_{bias} = a - cT^n$',fontsize=18)#  ,  $\frac{\chi^2}{DOF} = %.3f$'%red_chi_sq ,fontsize=18)
	pl.xlabel('Coldload Temperature [K]')
        pl.ylabel('Bias Power on TES [pW]')
	if legends:        
	    pl.legend()
	pl.grid(b=True)
	

