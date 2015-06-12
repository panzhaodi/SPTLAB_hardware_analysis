import pylab as pl
import numpy as np
import sptpol_software.util.files as files
import sptpol_software.util.time as time
from scipy import interpolate

rawfile='/home/cryo/Data/cryo/20141203_cooldown/log_2014-12-03T12h16.rawdat'
tempfile='/home/cryo/Data/cryo/20141203_cooldown/log_2014-12-03T12h16.dat'


resis = files.readColumns(rawfile, header = ['DATE', 'TIME', 'UCHEAD',  'ICHEAD',  'HEX',  'MAINPLATE',  'HE4PUMP', 'ICPUMP',
                                                  'UCPUMP',  'HE4SW',  'ICSW',  'UCSW',  'UC_STRAP',  'BOLO_CAN2',  'BOLO_CAN1',
						  'ICSTAGE',  '50K_PLATE',  '4K_PLATE',  'H_HE4SW', 'H_ICSW', 'H_UCSW',  
                                                  'H_HE4PUMP',  'H_ICPUMP',  'H_UCPUMP', 'H_COLDOAD','H_UCSTAGE', 'Frequency'], column_header = True)

temps = files.readColumns(tempfile, header = ['DATE', 'TIME', 'UCHEAD',  'ICHEAD',  'HEX',  'MAINPLATE',  'HE4PUMP', 'ICPUMP',
                                                  'UCPUMP',  'HE4SW',  'ICSW',  'UCSW',  'UC_STRAP',  'BOLO_CAN2',  'BOLO_CAN1',
						  'ICSTAGE',  '50K_PLATE',  '4K_PLATE',  'H_HE4SW', 'H_ICSW', 'H_UCSW',  
                                                  'H_HE4PUMP',  'H_ICPUMP',  'H_UCPUMP', 'H_COLDOAD','H_UCSTAGE', 'Frequency'], column_header = True)

date = [time.toSptDatetime(temps['DATE'][i].replace('/','-')+':'+temps['TIME'][i]) for i in range(0,len(temps['DATE']))] 

r=pl.array(resis['BOLO_CAN2'])
t=pl.array(temps['BOLO_CAN1'])

error=pl.where(t==-1)
r=delete(r,error)
t=delete(t,error)

r=list(r)
t=list(t)

r.sort(reverse=True)
t.sort()

interp_func=interpolate.interp1d(t,r)
print min(t),max(t)

file_temps=list(np.logspace(log10(min(t)),log10(max(t)),100))
print log10(min(t)),log10(max(t))
file_r=list(interp_func(file_temps))

file_temps.sort(reverse=True)
file_r.sort()

f=open('newout2.txt','w+')
f.write('Temp(K) R(Ohms)\n')

for i in range(len(file_temps)):
	string1 = str(file_temps[i])
	string2 = str(file_r[i])
	string = string1 + ' ' + string2 + '\n'
	f.write(string)
