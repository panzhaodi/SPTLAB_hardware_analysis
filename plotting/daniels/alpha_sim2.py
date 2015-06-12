from pylab import *
from math import erf

# Doesn't completely work

x=linspace(0,6,100)
figure(1)
clf()
Psat=9
Popt=linspace(0,8.99,1000)
r_fracs=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]

y=dict()

y['1']=arctan(x-3)-arctan(-3)
y['2']=tanh(x-3)-tanh(-3)
y['3']=(x-3)/(1+abs(x-3))+.75
y['4']=1/(1+e**(-x+3))-1/(1+e**(-3))
y['5']=[]
for i in x:
  y['5'].append(erf(i-3)-erf(-3))
y['5']=array(y['5'])

subplot(2,1,1)
plot(x,y['1']/y['1'][-1],label='arctan')
plot(x,y['2']/y['2'][-1],label='tanh')
plot(x,y['3']/y['3'][-1],label='abs')
plot(x,y['4']/y['4'][-1],label='log')
plot(x,y['5']/y['5'][-1],label='erf')
grid(b=True)
legend()

subplot(2,1,2)

alpha=dict()
alpha['1'] = (x*(y['1'][-1])/y['1'])/(x**2+1)  # T/R*dR/dT
alpha['2'] = (x*(y['2'][-1])/y['2'])*(-tanh(x)**2 + 1)  # T/R*dR/dT
alpha['3'] = (x*(y['3'][-1])/y['3'])*(-x*x/((abs(x) + 1)**2*abs(x)) + 1/(abs(x) + 1))  # T/R*dR/dT
alpha['4'] = (x*(y['4'][-1])/y['4'])*(e**x)/((e**x+1)**2)  # T/R*dR/dT
alpha['5'] = (x*(y['5'][-1])/y['5'])*(2*e**(-x**2))/sqrt(pi)  # T/R*dR/dT

plot(x,alpha['1'])
plot(x,alpha['2'])
plot(x,alpha['3'])
plot(x,alpha['4'])
plot(x,alpha['5'])
ylabel('Alpha')
xlabel('T')
grid(b=True)
figure(2)
clf()
plot(y['1']/y['1'][-1],alpha['1'])


for i in y.keys():  
  figure(3+(int(i)-1))
  clf()
  for j in range(len(r_fracs)):
      vbias2 = (Psat-4.5)*r_fracs[j]
      R = vbias2/(Psat-Popt)
      rf = R[where(R< 1.0)[0]]
      a=interp(rf,y[i]/y[i][-1],alpha[i])
      loopgain=a*(Psat-Popt[0:len(rf)])/(3*Psat)
      plot(Popt[0:len(rf)],loopgain,label='Rfrac='+str(r_fracs[j]))
      xlabel('Total Optical Power [pW]')
      ylabel('Loopgain')
      title(str(i)+' biased with Popt = 4.5 pW')
      grid(b=True)
      legend()

