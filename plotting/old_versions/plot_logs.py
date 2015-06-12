import pylab as pl
import numpy as np
import sptpol_software.util.files as files
import sptpol_software.util.time as time

def temp_plot(devices, others=['UCSTAGE','ICSTAGE'],log_file='/home/cryo/Data/cryo/20131121_cooldown/SCTemp.20131121.stoli1.log',
              one_plot=True, figure=None,ic_uc_temps='/home/cryo/Data/cryo/20140213_cooldown/IC_UC_temps.log'):
    '''
    Plot the temperature vs time of the devices you specify.                 

    device options are: 
        (if you give more than one, give them as a list)
           For Thermometry:
                'UCHEAD', 'ICHEAD', 
                'UCSTAGE', 'ICSTAGE', 
                'MAINPLATE', 'HEX',
                '4HEPUMP', 'ICPUMP', 'UCPUMP', 
                '4HESW', 'ICSW', 'UCSW', 
                'DIODEB1', 'DIODEB2', 'DIODEB3', 'DIODEB4'
           For Heater Currents:
                'H_4HEPUMP', 'H_ICPUMP', 'H_UCPUMP'
           For Boolean Switches:
                'S_4HESW', 'S_ICSW', 'S_UCSW'
    
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
    temps = files.readColumns(log_file, header = ['DATE', 'TIME',
                                                  'UCHEAD', 'ICHEAD', 
                                                  'UCSTAGE', 'ICSTAGE', 
                                                  'MAINPLATE', 'HEX',
                                                  '4HEPUMP', 'ICPUMP', 'UCPUMP', 
                                                  '4HESW', 'ICSW', 'UCSW', 
                                                  'DIODEB1', 'DIODEB2', 'DIODEB3', 'DIODEB4', 
                                                  'H_4HEPUMP', 'H_ICPUMP', 'H_UCPUMP', 
                                                  'S_4HESW', 'S_ICSW', 'S_UCSW'],
                              column_header = True)

    date1 = [time.toSptDatetime(temps['DATE'][i].replace('/','-')+':'+temps['TIME'][i]) for i in range(0,len(temps['DATE']))]
#    date = [time.datetime_to_mjd(temps['DATE'][i].replace('/','-')+':'+temps['TIME'][i]) for i in range(0,len(temps['DATE']))]

    #and the other log file
    if others:    
	ic_uc = files.readColumns(ic_uc_temps, header = ['DATE', 'TIME', 'UCSTAGE', 'ICSTAGE'],
                              column_header = True)
	date2 = [time.toSptDatetime(ic_uc['DATE'][i].replace('/','-')+':'+ic_uc['TIME'][i]) for i in range(0,len(ic_uc['DATE']))]

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

        pl.plot_date(date1,temps[_dev],fmt='-',label=_dev)
	#print min(temps['UCHEAD'][11982:])
        
    if others:
      for _dev in others:
	pl.plot_date(date2,ic_uc[_dev],fmt='-',label=_dev)
	pl.xlabel('Time')

    pl.legend()
    pl.grid()

#temp_plot(['UCHEAD'], others=False)


        
    
