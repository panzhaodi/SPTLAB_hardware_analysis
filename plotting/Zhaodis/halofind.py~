import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import fft as fft
import pickle
fname = '/home/cryo/Code/uofc/plotting/Zhaodis/lcdm64_step_step3000'

# read snapshot:
snapshot = open(fname + '.dat', 'r')

[Lbox,a] = np.fromfile(snapshot, np.float32, 2)
[step, Ngrid, Npart] = np.fromfile(snapshot, np.int32, 3)
x=np.fromfile(snapshot, np.float32, Npart)
y=np.fromfile(snapshot, np.float32, Npart)
z=np.fromfile(snapshot, np.float32, Npart)
px=np.fromfile(snapshot, np.float32, Npart)
py=np.fromfile(snapshot, np.float32, Npart)
pz=np.fromfile(snapshot, np.float32, Npart)
phi=np.fromfile(snapshot, np.float32, Npart)
rho=np.fromfile(snapshot, np.float32, Ngrid**3).reshape(Ngrid,Ngrid,Ngrid)

particles =[]
halos=[]
halo=[]
#set up list to store particle information
for i in range(0,Npart):
   particles.append([i,x[i],y[i],z[i], phi[i], 0.0])

#this is a loop for halos. j is the halo index
for j in range(0,632):
        length=len(particles)
        halo=[]
	imin=0; phimin=particles[0][4]
#find the minimum potential
	for i in range(0,length):
  	  if particles[i][4]<phimin:
   	      phimin=particles[i][4]; imin=i
        xcore=particles[imin][1]; ycore=particles[imin][2]; zcore=particles[imin][3]

#calculate the distance between every particle and this one
	for i in range(0,length):
             particles[i][5] = np.sqrt(np.power((particles[i][1]%Ngrid-xcore%Ngrid),2)+np.power((particles[i][2]%Ngrid-ycore%Ngrid),2)+np.power((particles[i][3]%Ngrid-zcore%Ngrid),2))


        particles=sorted(particles, key=lambda particles: particles[5])
        density = 10000
        for i in range(0,length):
             distance=particles[i][5]
             volume=4.0/3.0*np.pi*distance**3
             if not i==0:
                density = (1.0+i)/(volume)          
             if density < 731.34 or distance > 4.0/1.28 :
                break
             #store the list in halo[]
             halo.append(particles[i][0])
             #if i<10: print particles[i][0]
        print halo
        print j
        #delete list elements
        particles=particles[len(halo)+3000:]
        halos.append(halo)
with open('halo.pkl', 'wb') as f:
    pickle.dump(halos, f)


	

