from pylab import *

# Trying to see the effect of the shape of an RT curve on alpha.
# Alpha is directly proportional to the "slopes", so its a little redundant.

x=linspace(0,2.5*pi+8,1000)
figure(1)
Psat=9
Popt=linspace(0,8.99,1000)
slopes=[1,2,10,100]
r_fracs=[0.3,0.4,0.5,0.6,0.7,0.8,0.9]

clf()

figure(2)
clf()

figure(3)
clf()

for i in slopes:
  y=arctan(i*(x-8))-arctan(-i*8)
  figure(1)
  subplot(2,1,1)
  plot(x,y/y[-1],label = 'arctan('+str(i)+'x)')
  ylabel('Resistance')
  alpha = (x*y[-1]/y)*(i/(y[-1]*((8-x)**2)*(i**2) +1))  # T/R*dR/dT
  subplot(2,1,2)
  plot(x,alpha)
  ylabel('Alpha')
  xlabel('T')
  grid(b=True)
  figure(2)
  plot(y/y[-1],alpha)
  
  figure(3+slopes.index(i))
  clf()
  for j in range(len(r_fracs)):
      vbias2 = (Psat-4.5)*r_fracs[j]
      R = vbias2/(Psat-Popt)
      rf = R[where(R< 1.0)[0]]
      a=interp(rf,y/y[-1],alpha)
      loopgain=a*(Psat-Popt[0:len(rf)])/(3*Psat)
      plot(Popt[0:len(rf)],loopgain,label='Rfrac='+str(r_fracs[j]))
      xlabel('Total Optical Power [pW]')
      ylabel('Loopgain')
      title('Steepness level '+str(i)+' biased with Popt = 4.5 pW')
      grid(b=True)
      legend()
  
figure(1)
subplot(2,1,1)
legend(loc=4)
grid(b=True)

figure(2)
legend(loc=4)
grid(b=True)
ylabel('Alpha')

