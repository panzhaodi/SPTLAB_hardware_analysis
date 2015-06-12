import pylab as pl
import numpy as np
import cPickle as pkl


def plotEtf(item):
    
    try:
        item.keys()
    except AttributeError:
        print 'reading pickle'
        # if it does not have keys, assume it is a pkl file 
        #   and read it in
        item = pkl.load(open(item,'r'))

    chans = [i for i in item.keys() if i !='status' and i!='arguments']
    for chan in chans:
        freqs = item[chan]['frequencies']
        freqs_lab = [str(int(freq))+' Hz' for freq in freqs]
        
        sines = ([i for i in item[chan].keys() if i[-2:]=='Hz'])
        temp = zip([int(i[:-3]) for i in sines],sines)
        temp.sort()
        sines,junk = zip(*temp)

        ydatas = np.sort([i for i in item[chan].keys() if i[0:2]=='y1'])
        
        if len(ydatas)<len(freqs):
            print 'There is no data for frequencies: ' + str(freqs_lab[len(freqs):])
            freqs = freqs[:len(freqs)]
            
            
        pl.figure(20)
        pl.clf()
        for _freq,_freq_lab,_ydata_lab in zip(freqs,freqs_lab,ydatas):
            _ydata = item[chan][_ydata_lab]
            fft = abs(np.fft.fft(_ydata))
            print 'assumes FIR6!!!'
            f_freq = np.fft.fftfreq(len(_ydata),d=1./(25.E6/2**17))
            pl.subplot(121)
            pl.plot(f_freq,fft,label=_freq_lab,alpha=0.6)
            pl.subplot(122)
            pl.plot(_ydata,label=_freq_lab,alpha=0.6)
        pl.legend()
            
