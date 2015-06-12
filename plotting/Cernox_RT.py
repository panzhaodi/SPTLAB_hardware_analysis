import pylab as pl
from scipy.optimize import curve_fit

temps=[]
resis=[]
#with open('/home/cryo/Data/cryo/SPT3G_UCstage_20131215.dat') as source_file:
with open('/home/cryo/Data/cryo/SPT3G_brick_20131215.dat') as source_file:
    for line in source_file.readlines()[3:]:
        cols = [float(x) for x in line.split()]
        temps.append(cols[0])
	resis.append(cols[1])
temps=pl.array(temps)[0:8]
resis=pl.array(resis)[0:8]*10**(-3)

#def fit_function(x,a,b,c,d,e):
	#return (a*x**4 + b*x**3 + c*x**2 + d*x + e)
def fit_function(x,a,b,c):
	return(a/(x**b) + c)

parameters=curve_fit(fit_function,resis,temps)
#[a,b,c,d,e]=parameters[0]
[a,b,c]=parameters[0]

pl.figure()
pl.scatter(resis,temps)
pl.ylabel('Temp [K]')
pl.xlabel('Resistance [kOhm]')
pl.grid()

x = pl.linspace(5,55,100)
#fit = fit_function(x,a,b,c,d,e)
fit=fit_function(x,a,b,c)
pl.plot(x,fit)

#print "[a,b,c,d,e]= ",[a,b,c,d,e]
print "[a,b,c]= ",[a,b,c]
pl.show(block=False)
