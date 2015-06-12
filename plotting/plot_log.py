# This plot_logs is currently somewhat inconvenient, as one must change it to match the header in the log file, which itself changes whenever the configuration changes. DPD 4/16/14

import pylab as pl
import numpy as np
import sptpol_software.util.files as files
import sptpol_software.util.time as time
from crython.cryoboard import CryoBoard


def temp_plot(devices, log_file='/home/cryo/Data/cryo/20150606_cooldown/log_2015-06-06T21h57.dat', one_plot=True, figure=None,return_data=False):
    '''
    Plot the temperature vs time of the devices you specify.                 

    device options are: 
        (if you give more than one, give them as a list)
           For Thermometry:
                'UCHEAD', 'ICHEAD', 
                'UCSTAGE', 'ICSTAGE', 
                'HEX', 'MAINPLATE',
                'HE4PUMP', 'ICPUMP', 'UCPUMP', 
                'HE4SW', 'ICSW', 'UCSW', 
                'UC_STRAP', 'BOLO_CAN2', 'BOLO_CAN1', 'IC_STAGE'
		'50K_PLATE','4K_PLATE'
               
    one_plot plots all given devices on the same plot

    figure is the number of the figure you want to plot on, 
            if it is the number of an open figure it will 
            clear the figure. If None then it will open the 
            next availible figure. 

    '''

    #Turn the devices into a list if it is not already
    if not isinstance(devices,list):
        devices = [devices]
    
    # read in the log file 
    temps = files.readColumns(log_file, header = ['DATE', 'TIME', 'UCHEAD',  'ICHEAD',  'HEX',  'MAINPLATE',  'HE4PUMP', 'ICPUMP',
                                                  'UCPUMP',  'HE4SW',  'ICSW',  'UCSW',  'WAFER_PLATE',  'IC_STAGE1',  'TC_BOX',
						  'IC_STAGE2',  'CL_TOP',  'CL_KEATING',  'H_HE4SW', 'H_ICSW', 'H_UCSW',  
                                                  'H_HE4PUMP',  'H_ICPUMP',  'H_UCPUMP', 'H_COLDOAD','H_UCSTAGE', 'Frequency'], column_header = True)

    date = [time.toSptDatetime(temps['DATE'][i].replace('/','-')+':'+temps['TIME'][i]) for i in range(0,len(temps['DATE']))] 

    if return_data:
	return date,temps  

    for _dnum,_dev in enumerate(devices):
        
        if _dnum==0:
            if figure:
                pl.figure(figure)
                pl.clf()
            else:
                pl.figure()
        else:
            if one_plot:
               pass
            else:
                pl.figure()

        pl.plot_date(date,temps[_dev],fmt='-',label=_dev)
	pl.xlabel('Time')
	pl.ylabel('Temp [K]')
    pl.legend()
    pl.grid()

        
    


