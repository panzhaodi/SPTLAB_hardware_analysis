import pylab as pl
import numpy as np
import cPickle as pkl
from scipy.interpolate import UnivariateSpline

# The oldest of the rt_plot family. Mainly written/used by Zhaodi, with first attempt at
# fitting alpha written by Daniel D.

def rt_plot(sq='2',plot_alpha=True):
  if sq =='2':
    document = '/home/cryo/Data/output/20140723_cooldown/tc139(copy)/20140731121808_139-1-0.dat'
    conversion_factor = 46.208
  if sq=='4':
    document = '/home/cryo/Data/output/20140723_cooldown/tc139(copy)/20140731121808_139-3-0.dat'
    conversion_factor = 45.596
  if sq=='5':
    document = '/home/cryo/Data/output/20140723_cooldown/tc137(copy)/20140731121946_137-0-0.dat'
    conversion_factor = 44.615
  if sq=='7':
    document = '/home/cryo/Data/output/20140723_cooldown/tc137(copy)/20140731121946_137-2-0.dat'
    conversion_factor = 45.278
  if sq =='8':
    document = '/home/cryo/Data/output/20140723_cooldown/tc137(copy)/20140731121946_137-3-0.dat'
    conversion_factor=  45.364
  with open(document) as f:
 	   data=f.read()
  data = data.split('\n')
  line = len(data)
  del data[0]
  del data[int(line-2)]
  line = len(data)
  x = [row.split('\t')[3] for row in data]
  y = [row.split('\t')[1] for row in data]
  for n in range (0, line):
    x[n]= 1000*float(x[n])
  for n in range (0, line):
    y[n]= conversion_factor/float(y[n])
    
  #colors=['.','.','r','.','g','b','.','k','m']
  pl.figure(1)
  pl.clf()
  pl.scatter(x,y,s=2)#,color=colors[int(sq)],label='Sq'+sq)
  pl.xlabel('Temperature (mK)')
  pl.ylabel('Resistence (Ohm)')
  pl.title('Sq'+sq+' RT measurement')
  pl.grid(b=True)

######### Fitting alpha #########
  if plot_alpha:
	    
    x=pl.sort(x)
    y=pl.sort(y)

    ##### Use only the data within the transition #####
    ind=pl.where((x > 514) & (x < 546))
    x=x[ind]
    y=y[ind] 

    pl.figure(4)
    pl.clf()
    pl.scatter(pl.log(x),pl.log(y),s=2)

    fit = UnivariateSpline(pl.log(x),pl.log(y),s=0.05,k=1)
    xfit = pl.linspace (6.245,6.30,1000)
    yfit=fit(xfit)
    
    pl.plot(xfit,yfit,'--')

    pl.figure(5)
    pl.clf()
    pl.subplot(2,1,1)
    pl.scatter(x,pl.array(y)/pl.sort(y)[-1],s=3)#,color=colors[int(sq)],label='Sq'+sq)
    pl.plot(pl.e**pl.array(xfit),pl.e**pl.array(yfit)/pl.e**pl.array(yfit)[-1],'-')
    pl.xlim(515,545)
    pl.ylabel('R/Rn')
    pl.suptitle('Sq'+sq+' RT measurement')
    pl.grid(b=True)

    alpha=[]
    xfittemp=[]
    for _x in xfit:
      alpha.append(fit.derivatives(_x)[1])
      xfittemp.append(pl.e**_x)
    pl.subplot(2,1,2)
    pl.plot(xfittemp,alpha,'-')
    pl.xlabel('Temperature (mK)')
    pl.ylabel('Alpha')
    pl.grid(b=True)


