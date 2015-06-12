import pylab as pl
from scipy.optimize import curve_fit
from collections import deque

# Written in IDL by L. Bleem, translated to Python by D. Dutcher. October 2014
# Contains all functions known to be working properly.


def where_closest(value,array):
  abs_diff = pl.array(abs(array-value))
  wh=where(abs_diff == min(abs_diff))[0]
  wh_closest = abs_diff[wh]

  return wh_closest




def create_band(band, res, plt=False, high=0, conv=False):
  '''
  print, "create_band, band, res, /bw, /plt, high=high, /conv"
  print, "This program will create a top hat band convolved with a gaussain of sigma = res with freq data spaced from 0 to 1000 (default) GHz"
  print, "band = bandpass region in GHz"
  print, "res = freq. spacing in GHz"
  print, "/plt plot bandpass"
  print, "high=high upper GHz region"
  print, "r = create_band([140, 160], 2.0, /bw, high=500.0)"
  return, 0
  '''

  npts = pl.ceil(4000.0)
  if high : npts = pl.ceil(high*1.0/.25)
  freq = pl.arange(npts)*.25 
  response = pl.zeros(len(freq))

  inb = pl.where((freq < band[1]) & (freq > band[0]))[0]
  if band[0] == band[1] : inb = pl.where_closest(freq, band(0))[0]
  if inb[0] == -1:
    print "Band not between 0-1000 GHZ"
    return 0

  response[inb] = 1

  #let's convolve the band with our resolution. 
  if conv:
    xx = .25*pl.arange(6*res/.25+1)-3*res
    con =pl.exp(-xx**2/(2*res**2))/pl.sqrt(2*pl.pi)/res
    normalization=1./sum(abs(con))
    response = pl.convolve(response, con,'same')*normalization
  

  if plt:
    pl.figure()
    pl.plot(freq, response,'D')
    pl.xlabel('Freq(GHz)')
    pl.ylabel('Reponse')
    pl.xlim(band[0] - 3*res, band[1] + 3*res)


  result = {'Freq': freq, 'resp': response}

  return result




def create_interf(freq,resp,band=[], plt=False,sav=False, res=1.0, two=False):
  '''
  print, "create_interf, freq, resp, tc=tc, plt=plt,sav=sav, band=band, res=res, bw=bw, two=two"
  print, "freq, resp - put in your own frequency and response data"
  print, "/plt plots the band pass and interferrogram"
  print, "/sav saves the interferrogram to a text file"
  print, "band = band, res=res, /bw - use these to create freq/resp band with create_band"
  print, "/two - put 2 interferrograms in a row"
  return, 0
  '''


#  def where_closest(value,array):
#    abs_diff = pl.array(abs(array-value))
#    wh=pl.where(abs_diff == min(abs_diff))[0]
#    wh_closest = abs_diff[wh]
#    return wh_closest

  if len(band) != 0:
    r = create_band(band, res)
    freq = r['Freq']
    resp = r['resp']
    if band[1] == band[0]:
      resp = pl.zeros(len(freq))
      k = where_closest(band[0], freq)
      resp[k[0]] = 1.0

#if freq(0) != 0 then return with warning!
  if freq[0] != 0:
    print 'Must go down to zero frequency'
    return -1

  #Let's be careful with these n's
  #From DFanning
  #Let N=8, N/2 = 4, then F_ns for 0,1,2,3,4, -3,-2,-1  NOTE: no -4!

  n = pl.arange(len(freq)/2.+1)
  x = 30*n/(max(freq)-min(freq)) #30 to go from GHz to icm


  intf = pl.ifft(resp)
  x2 = pl.concatenate((x, -(x[1:len(x)-2])[::-1]))    #Crap. should this be -2 or -1
  if len(freq) % 2 == 1 : x2 = pl.concatenate((x, -(x[1:len(x)-1])[::-1])) 

  #plot, freq, resp
  #oplot, freq, FFT(intf), color=1

  if two:
    x2 = pl.concatenate((x2, x2+2*max(x2)))
    intf = pl.concatenate((intf, intf))

  q = x2.argsort()
  x2 = x2[q]
  intf = (intf[q]).real
  result ={'x': x2, 'intf':intf}

  if plt:
    if len(band) != 0 : rtemp = create_band(band, res, plt=True)
    pl.plot(freq, resp)
    pl.title('Band')
    pl.figure()
    pl.plot(x2, intf.real)
    pl.xlabel('Optical Delay (cm)')
    pl.ylabel('Response') 
    pl.title('Interferrogram')


#if sav:
#   openw, 1, sav
#   x0 = result.x(0:n_elements(x2)-1)
#   outp = [[x0], [real_part(result.intf)], [imaginary(result.intf)]]
#   printf, 1, transpose(outp)
#   close, 1

  return result




def find_intf(cts, intf, plt=False):

  intf2 = intf/max(intf)

  width = len(intf)
  conv = pl.zeros(len(cts) - width)

  #Put this in histogram form instead of for-loop form!
  for i in range(0,len(cts)-width):
    xx = cts[i:i+width] -(cts[i:i+width].mean())
    conv[i] = (sum((intf2*xx)))

  spacer = pl.zeros(pl.floor(width/2.0))
  conv = pl.concatenate((spacer, conv))


  output = {'conv': conv, 'index': pl.arange(len(conv)), 'shift': pl.floor(width/2.0)}


  if plt:
    pl.figure()
    pl.plot(output['index'], output['conv'], 'D') 
    pl.title('Convolution vs. index')
    pl.xlim(0, len(cts))
    pl.xlabel('Index')
    pl.figure()
    pl.plot(range(len(cts)),cts)
    pl.xlim(0, len (cts))
    pl.title('Time stream')


  return output





def find_fts_peaks_max(conv, optv, rate, thres, length, adc_cts2, onepk=False, nscans=1e6):

  number_intf = 2*nscans

  maxpk = max(conv['conv'])
  pos = pl.arange(len(conv['conv']))*optv/rate 
  iid = pl.arange(len(conv['conv'])) 
  conv2 = conv['conv']
  adc_cts = adc_cts2
  xax = pl.arange(len(adc_cts))
  space = .2     #in cm to look around conv. peak in the actual time stream data
               #for the "true" white light fringe
  gpk = 0

  #First find the first max of the data
  mm = pl.where(conv2 == max(conv2))[0]
  index = iid[mm[0]]
  #Look within "space" around this peak for the maximum of ADC cts
  small_reg = pl.where((pos <= (pos[mm[0]] + space)) & (pos >= (pos[mm[0]] - space)))[0]
  adc_pk = pl.where(adc_cts[small_reg] == max(adc_cts[small_reg]))[0]
  index = iid[small_reg[adc_pk[0]]]
  cindex = [iid[mm[0]]]

  if onepk: #only looking for 1 peak, so return this first one
    output = {'tpks': index, 'cpks': cindex}
    return output

  #we will looks for peaks with outliers in adc cts. Need to deal with the slope.
  rt = pl.polyfit(xax, adc_cts,1)   #first try - removing a linear fix from the whole scan
  adc_cts = adc_cts - rt[1] - rt[0]*xax
  mm = [small_reg[adc_pk[0]]]
  madc = adc_cts[mm[0]]

  while gpk < 1 :    #putting loop in to try and protect against DC spikes, only positive spikes mess up the peak finder
  #Next cut out in length around this pt
    newr = pl.where((pos < (pos[mm[0]] - length)) | (pos > (pos[mm[0]] + length)))[0]
    conv2 = conv2[newr]
    pos = pos[newr]
    iid = iid[newr]
    adc_cts = adc_cts[newr]

    #again find the peak in adc cts:
    mm = pl.where(conv2 == max(conv2))[0]
    small_reg = pl.where((pos <= (pos[mm[0]] + space)) & (pos >= (pos[mm[0]] - space)))[0]
    adc_pk = pl.where(adc_cts[small_reg] == max(adc_cts[small_reg]))[0]
    tadc_new = adc_cts[small_reg[adc_pk[0]]]

#print, conv2(mm(0)), .75*conv.conv(cindex(0))
    initial = .75*conv['conv'][cindex[0]]
    testcase = conv2[mm[0]]
    if testcase > initial:  #The first peak is Good!
       gpk = 1.0
       cindex = pl.concatenate((cindex, [iid[mm[0]]]))
       index = pl.array([index, iid[small_reg[adc_pk[0]]]])
       newr = pl.where((pos < (pos[mm[0]] - length)) | (pos > (pos[mm[0]] + length)))[0]

       if len(newr) == 0:
         #check that we are above our threshold cutoff
         accept_pks = pl.where(conv['conv'][cindex] >= thres*conv['conv'][cindex[0]])[0]
         output = {'tpks': index[accept_pks], 'cpks': cindex[accept_pks]}
         return output
       conv2 = conv2[newr]
       pos = pos[newr]
       iid = iid[newr]
       adc_cts = adc_cts[newr]
   
    if testcase < initial: #It's Bad!
      madc = tadc_new
      maxpk = max(conv2)
      cindex = iid[mm[0]]
      index = iid[small_reg[adc_pk[0]]]


  pk = max(conv2)

  while ((pk >= (thres*maxpk)) & (len(index) < number_intf)) :
    mm = pl.where(conv2 == max(conv2))[0]
    cindex = pl.concatenate((cindex, [iid[mm[0]]]))
    pk = conv2[mm[0]]
    small_reg = pl.where((pos <= (pos[mm[0]] + space)) & (pos >= (pos[mm[0]] - space)))[0]
    adc_pk = pl.where(adc_cts[small_reg] == max(adc_cts[small_reg]))[0]
    mm = [small_reg[adc_pk[0]]]
    index = pl.concatenate((index, [iid[small_reg[adc_pk[0]]]]))
    newr = pl.where((pos < (pos[mm[0]] - length)) |  (pos > (pos[mm[0]] + length)))[0]
    if len(newr) != 0:
      conv2 = conv2[newr]
      pos = pos[newr]
      iid = iid[newr]
      adc_cts = adc_cts[newr] 
    else: pk = 0   

  pts = len(index)
  #check that we are above our threshold cutoff
  accept_pks = pl.where(conv['conv'][cindex] >= thres*conv['conv'][cindex[0]])[0]
  output = {'tpks': index[accept_pks], 'cpks': cindex[accept_pks]}


  return output




def fts_fft(cts, optv,rate, pks2, length, band=[], hann=False, bolo=0, plt=False,absv=False, phaseb=[], chop=False, crange=[], notquiet=False):
  '''
  IF N_PARAMS() == 0:
   print('pro fts_fft, cts, optv,rate, pks2, length, band=band, hann=hann, bolo=bolo, abs=abs, phaseb=phaseb, chop=chop,   crange=crange, notquiet=notquiet'
   print('To be added TC deconvolution'
   print('/abs -> phase corrected interferrogram is from both real and imaginary FFT otherwise phase corrected interferrogram just the real part.'
  ENDIF
  '''
  #Document the output
  #result = {'freq': freq_a, 'real_FFT': real_a, 'im_FFT': im_a,'abs': abs_a, 'time_scan': timeout, 'whitel': pks, 'xint':   xint_a, 'intf': int_a, 'scan_length':scanl}

  if len(phaseb)==0: phaseb=band

  pks = pl.sort(pks2)
  pos = pl.arange(len(cts))*optv/rate
  time = pks/rate

  for i in range(len(pks)):

                                #grab the data within length of each
                                #white light peak and make scan
                                #symmetric around peak
    dd = pl.where((pos > (pos[pks[i]] - length)) & (pos<pos[pks[i]] + length))[0]
    d1 = pl.where(pos[dd] > (pos[pks[i]]))[0]
    d2 = pl.where(pos[dd]<(pos[pks[i]]))[0]

    mmin = min([len(d1), len(d2)]) #figuring out if one half is shorter than the other
    cts_range = pl.arange(2*mmin) + pks[i]-mmin
    #take out the mean and slope:
    xxtemp = pl.arange(len(cts_range))
    rrtemp = pl.polyfit(xxtemp, cts[cts_range],1)
    yr = cts[cts_range] - rrtemp[0]*xxtemp
    yr = yr - yr.mean()
    #Need even length for FFT
    if len(yr) % 2 != 0 : yr = yr[1:len(yr)]

#deal with possible chopping now:
#basic idea: find the chopper peak frequency - then we'll take
#                                              the FFT, and set this
#chopper peak to 0 frequency. Use the negative side band as our signal (lose 1/2
#the signal, but then don't have to deal with potentially large
#change in frequency response between the two side bands). 


    if chop != 0:

      if len(crange) == 0:
        print(crange)
        chophi = 15
        choplow = 8
      else:
        chophi = crange[1]
        choplow = crange[0]
  
 #fitting chopped signal FFT with a gaussian around the chopper peak 

      def gaussian(x,a0,a1,a2,a3): # Mimics IDL's gaussfit with 'nterms' = 4
        z = (x-a1)/a2
        y = a0*pl.exp(-(z**2)/2) + a3
        return y
        

      qout = time_fft(yr, samplerate =rate, hann=True)                              ######Does this work?#######
      q2 = time_fft((qout['real']+1j*qout['im']), inverse=True)
      fr = pl.where((qout['freq'] > choplow) & (qout['freq']<chophi))[0]
      fit, pcov = curve_fit(gaussian,qout['freq'][fr],qout['abs'][fr])       ## This fit has inf covariance matrix. Not great.
      pkf = fit[1]  #found the peak response frequency.
      chspec_f = 30.0*pkf/optv  #this is where it maps in GHz
      if notquiet:
        print("Chopper at "+str(pkf)+ " Hz")
        print("Spectra: "+str(chspec_f)+ " GHz")
        print('3rd Harmonic in subtracted scan at ',+str(2*chspec_f)+" GHz")

#Save orignal for comparision with created intf
    yout = yr

    if hann:
      w1 = pl.hanning(len(yr))
      yr = yr*w1

    #Let's get this shifting right:
    yr=deque(yr)
    yr.rotate(int(-len(yr)/2.0 +1))
    yr=pl.array(yr)
    n = pl.arange(len(yr)/2. + 1)
    n2 = pl.concatenate((n, -(n[1:len(n)-1])[::-1]))    #CRAP should this be -2 or-1?
    icm = n2/len(yr)/(optv/rate)
    icm0 = icm
    FFT_r = pl.fft(yr)

    icm2 = icm
    #okay now let's do the shifting and such with the chopper:
    if chop != 0:
      chspec_icm = chspec_f/30.0
      icm = -1.0*(icm-chspec_icm)

    if notquiet:
      pl.figure()
    
      pl.plot(30*icm, abs(FFT(yr)), label='Not Demodulated')
      pl.plot(30*icm2, abs(FFT(yr)), label='Lower Side Band')
      pl.plot(-30*icm, abs(FFT(yr)), label='Upper Side Band')
      pl.xlim(50, 300) 
      pl.title('Chopped data, abs value')
      pl.ylim(0, .5)
   #stop


   #now to take out a phase:
    phase = pl.arctan2(FFT_r.imag, FFT_r.real)
    tt =pl.where((icm > phaseb[0]/30.) & (icm<phaseb[1]/30.))[0]

    def linfit(x,a,b):
      y=a+b*x
      return y    

    if len(tt) == 1 : r = [0,0] 
    if len(tt) != 1 : r,pcov = curve_fit(linfit,icm[tt], phase[tt],sigma = 1/abs(FFT_r[tt])**2)  # 'sigma' here was 'measure_errors' in IDL

    #apply phase correction:

    pc = r[0] + r[1]*icm
    shift_c = pl.exp(-pc*1j)
     
    FFT_r = FFT_r*shift_c
    FFT_net = FFT_r


    #create the interferrogram to feed back:
    intf = create_interf(30*icm0,pl.fft(yr).real)
    if absv : intf = create_interf(30*icm, pl.fft(yr))
    xin = intf['x']
    intfa = intf['intf']
    length_int =[len(intfa)]
    #pl.plot(intf.x, intf.intf, xr = [-2, 2]


    #keep only positive frequencies:
    qq = pl.where(icm >= 0)
    icm = icm[qq]
    FFT_r = FFT_r[qq]

    #sort this
    qqsort = icm.argsort()
    icm = icm[qqsort]
    FFT_r = FFT_r[qqsort]


    if len(band) == 0 : band = [50, 700]

    ptitle = 'FTS Data for bolo ' + str(bolo)
    if i == 0:
      freqout = 30.0*icm
      realout = (FFT_r[0:len(icm)]).real
      imout = (FFT_r[0:len(icm)]).imag
      sindex = pl.zeros(len(icm))+i
      timeout = [time[i]]
      length_inter = length_int
      xint = xin
      intfg = intfa
      #stop     
      if plt:
         pl.figure()
         inband = pl.where((30*icm > band[0]) & (30*icm < band[1]))[0]
         if len(inband) == 0 : inband = pl.arange(len(icm))
         if absv == 0:
           pl.plot(30*icm[1:], (FFT_r[1:]).real/max((FFT_r[inband]).real),'k-',label='%.4f cm'%(time[i]*optv))
           pl.plot(30*icm[1:], (FFT_r[1:]).imag/max((FFT_r[inband]).real), 'k--')
           pl.xlabel('Frequency (GHz)')
           pl.ylabel('Normalized Spectra')
           pl.xlim(band[0],band[1])
           pl.title(ptitle)
           pl.ylim(-0.5, 1)

         #stop
         if absv != 0:
              pl.plot(30*icm[1:], abs(FFT_r[1:])/max(abs(FFT_r[inband])))
              pl.xlabel('Frequency (GHz)')
              pl.ylabel('Normalized Abs Value of Spectra')
              pl.xlim(band)
              pl.title(ptitle)
              pl.xlim(-.1, 1)

    if i != 0:
      tfreqout = 30*icm
      trealout = (FFT_r[0:len(icm)]).real
      timout = (FFT_r[0:len(icm)]).imag
      tsindex = pl.zeros(len(icm))+i      
      ttimeout = [time[i]]

      freqout = pl.concatenate((freqout, tfreqout))  ## Don't know if these are lists,arrays or integers
      realout = pl.concatenate((realout, trealout))
      imout = pl.concatenate((imout, timout))
      sindex = pl.concatenate((sindex, tsindex))
      timeout = pl.concatenate((timeout, ttimeout))
      length_inter = pl.concatenate((length_inter, length_int))
      xint = pl.concatenate((xint, xin))
      intfg = pl.concatenate((intfg, intfa))


      if plt:
        inband = pl.where((30*icm > band[0]) & (30*icm<band[1]))[0]
        if inband[0] == -1 : inband = pl.arange(len(icm))
        if absv == 0:
          pl.plot(30*icm[1:], (FFT_r[1:]).real/max((FFT_r[inband]).real),color=pl.cm.jet(.15*i),label='%.4f cm'%(time[i]*optv))
          pl.plot(30*icm[1:], (FFT_r[1:]).imag/max((FFT_r[inband]).real), '--',color=pl.cm.jet(.15*i))
        if absv != 0:
          pl.plot(30*icm[1:], abs(FFT_r[1:])/max(abs(FFT_r[inband])),color=pl.cm.jet(.15*i),label='%.4f cm'%(time[i]*optv))
        pl.legend(loc=4)

  pl.grid()


#okay, let's get this result in some resasonable form:

  scanl = pl.histogram(sindex, bins = int(max(sindex)-min(sindex)+1))[0] #figure out the lengths of each scan and pad array with zeros if necessary
  maxl = max(scanl)
  scan_index = pl.arange(len(pks))
  freq_a = pl.zeros((maxl, len(pks)))           ##IDL is backwards, had to flip all these from CxR to RxC
  real_a = pl.zeros((maxl, len(pks)))
  im_a = pl.zeros((maxl, len(pks)))
  abs_a = pl.zeros((maxl, len(pks)))

  for i in range(len(pks)):  ## The below all comes from the column-focused IDL. Re-write? ##

    if i == 0:
      if maxl - scanl[i] != 0:  
        zeros_to_add = pl.zeros(maxl - scanl[i])             ## Don't know if these are arrays or integers ##
        freq_a[:,i] = pl.concatenate((freqout[0:scanl[i]],zeros_to_add))
        real_a[:,i] = pl.concatenate((realout[0:scanl[i]],zeros_to_add))
        im_a[:,i] = pl.concatenate((imout[0:scanl[i]],zeros_to_add))
        abs_a[:,i] = pl.sqrt(real_a[:,i]**2 + im_a[:,i]**2)
        new_start = scanl[i]
      else: 
        freq_a[:,i] = freqout[0:scanl[i]]
        real_a[:,i] = realout[0:scanl[i]]
        im_a[:,i] =imout[0:scanl[i]]
        abs_a[:,i] = pl.sqrt(real_a[:,i]**2 + im_a[:,i]**2)
        new_start = scanl[i]
      
    else:
      if maxl - scanl[i] != 0:  
        zeros_to_add = pl.zeros(maxl - scanl[i])
        freq_a[:,i] = pl.concatenate((freqout[new_start:new_start +scanl[i]],zeros_to_add))
        real_a[:,i] = pl.concatenate((realout[new_start:new_start +scanl[i]],zeros_to_add))
        im_a[:,i] =pl.concatenate((imout[new_start:new_start + scanl[i]],zeros_to_add))
        abs_a[:,i] = pl.sqrt(real_a[:,i]**2 + im_a[:,i]**2)
        new_start = new_start + scanl[i]
      else: 
        freq_a[:,i] = freqout[new_start:new_start +scanl[i]]
        real_a[:,i] = realout[new_start:new_start +scanl[i]]
        im_a[:,i] =imout[new_start:new_start +scanl[i]]
        abs_a[:,i] = pl.sqrt(real_a[:,i]**2 + im_a[:,i]**2)
        new_start = new_start + scanl[i]

#now to deal with the interferrograms:

  maxl = max(length_inter)
  xint_a = pl.zeros((maxl, len(length_inter)))
  int_a = pl.zeros((maxl, len(length_inter)))

  for i in range(len(pks)):

    if i == 0:
      if maxl - length_inter[i] != 0:
        zeros_to_add = pl.zeros(maxl - length_inter[i])
        xint_a[:,i] = pl.concatenate((zeros_to_add, xint[0:length_inter[i]]))
        int_a[:,i] = pl.concatenate((zeros_to_add, intfg[0:length_inter[i]]))
        new_start = length_inter[i]
      else:
        xint_a[:,i] = xint[0:length_inter[i]]
        int_a[:,i] = intfg[0:length_inter[i]]
        new_start = length_inter[i]

    else:
      if maxl - length_inter[i] != 0:
        zeros_to_add = pl.zeros(maxl - length_inter[i])
        xint_a[:,i] = pl.concatenate((zeros_to_add, xint[new_start:new_start+length_inter[i]]))
        int_a[:,i] = pl.concatenate((zeros_to_add, intfg[new_start:new_start+length_inter[i]]))
        new_start = new_start +  length_inter[i]
      else:
        xint_a[:,i] = xint[new_start:new_start+length_inter[i]]
        int_a[:,i] = intfg[new_start:new_start+length_inter[i]]


  result = {'freq': freq_a, 'real_FFT': real_a, 'im_FFT': im_a,'abs': abs_a, 'time_scan': timeout, 'whitel': pks, 'xint': xint_a, 'intf': int_a, 'scan_length':scanl}



  return result







def analyze_fts_scan( adc_cts, optv, band=[145, 178], rate=191, thres=.5, length=10, sav=False, bolo='', hann=False, pband=[], plot=True, exclude=[], include=[],  absv=False,template=False, phaseb=[],chop=False, chrange = [14,16], onepk=False, nscans=1e6, notquiet=False):


#IF N_PARAMS() == 0:
#   print, "outp= analyze_fts_scan(adc_cts, optv, band=band, rate=rate, thres=thres, length=length, sav=sav, bolo=bolo, hann=hann, pband=pband, noplt=noplt, #exclude=exclude, include=include,  abs=abs,template=template, batch=batch, phaseb=phaseb, back=back, chop=chop, chrange = chrange, onepk=onepk, nscans=nscans, #notquiet=notquiet)"
#   return
#

  '''
  #+
  # NAME:
  #   ANALYZE_FTS_SCAN
  # PURPOSE:
  #       Find and analyze interferograms from a timestream
  # CALLING SEQUENCE:
  #       output = analyze_fts_scan(data, optv, band=band, bw=bw, rate=rate,
  #       res=res, thres=thres, length=length, sav=sav,
  #       bolo=bolo, hann=hann, pband=pband, plt=plt, exclude=exclude, include=include)
  # EXAMPLE:
  #       outp = analyze_scan(timestream, 2.0, bolo=415)
  #       outp = analyze_scan(timestream, 2.0, band=[195,220],
  #       rate=1000., length=13.0, bolo="2a:ra2", /hann, exclude = [1,
  #       2], chop=1, chop_range = [14, 16])
  # KEYWORDS:
  #       data - the timestream data (1D array)
  #       optv - optical velocity
  # OPTIONAL KEYWORDS:
  #       template - interferrogram template, IDL save file with
  #                  variables tempx, tempy, xranges from -.8 to .8 cm
  #       band - [Start GHz, Stop GHz], the band of the template
  #              interferrogram used to find the white light peaks in
  #              the time stream. By default the template is a flat band
  #              for 150GHz
  #       rate - sampling rate of the data, default is 100Hz  
  #       thres - optional,fraction of max convolution peak to label as
  #               a white light fringe, mess with this if you are not
  #               finding all the white light peaks in the data, or if
  #               the amplitudes of the interferograms drastically
  #               change mid time stream
  #       length - maximum optical delay to include when taking FFT of
  #                interferogram, default 10 cm
  #       /sav - writes out an idl save file, bolo_N_processed_FFT
  #       bolo - put in name/number/description of bolometer/scan in
  #              this parameter, used in save file name, naming plot
  #              windows
  #       /hann - apply a hanning window to data before FFT
  #       pband - [Start GHz, Stop GHz] band to plot, default 50-800 GHz
  #       /plt - If you want to plot the FFT data as well as the
  #              timestream/white light locator (default no plots at all)
  #      /abs - Plot the absolute values of the interferrograms, use for
  #             really noisy data     
  #      exclude - array containing the numbers of scans to exclude
  #                 from plots (scan numbering starts at 1)
  #      include - array containing the numbers of scans to include,
  #                 scan numbering starts at 1
  #      /noise - take the maximums of the correlation functions as the
  #              fit of a gaussian around the peak - commented out in code
  #      /batch - call this when you are batch processing a bunch of FFT's
  #      phaseb = [start GHz, stopGHz], region in which to fit out a
  #      linear phase shift (by default uses the entire FFT region with
  #      error bars (1/abs(FFT)^2)
  #      chop = set this to a non zero number if the interferrogram is
  #      riding on a chopped signal, by default chop must be between 8
  #      and 15 Hz, other wise set chop_range. NOTE: INFORMATION IN YOUR
  #      BAND SHOULD BE GREATER THAN .5 HZ AWAY FROM THE CHOP SIGNAL!!
  #      Otherwise, demodulation as currently done (12/2/08) will mess up your signal
  #      chrange = [lower, upper], range of frequencies in which to
  #      look for the chopper signal
  #      onepk=onepk   set this not equal to zero if there is only 1
  #      interferrogram in the data
  #      /thick -> include if you are chopping with the thick grill
  #      nscans = # of scans (2 interferrograms/scan)
  #
  # OUTPUT:
  #       output a structure with the following entries:
  #       Freq - freq data, the scans are one after another
  #       real_fft - real part of scans FFT
  #       im_fft - imaginary " "
  #       abs - absolute value of FFT  

  #       time_scan - time white light peak occured after start of file
  #                   - if these are oddly spaced, prob. missed an
  #                     interferogram or found a spurious one
  #       to plot a particular scan, 
  #          plot, output.freq(scan, *), output.real_part(scan, *), xr = [50, 700]
  #       whitel - constains indexes of found white light peaks in time stream
  #       xint - cm for phase corrected interferrogram
  #       intf - phase corrected interferrogram. Is the FFT of the real
  #              part of the phase corrected Spectrum. HOWEVER, if call
  #              /abs then this is the FFT of both the real and
  #              imaginary part of the phase corrected interferrogram. 
  # TO DO:
  #     Right now, only removes a linear phase term in the FFT, would be
  #     nice to add a time constant fitter or deconvolver.
  # CREATED:
  #     Nov. 19, 2008: LB
  # Converted to python:
  #     2014 DPD
  #-
  '''


  #calls find_intf, create_intf, create_band, find_fts_peaks, fts_FFT
  chop1=chop
  adc_cts2 = adc_cts

  chophi = chrange[1]
  choplow = chrange[0]


#-----------------------------------------
#What to do about the chopped signal? Right now just using spurious
#reflections to find the interferrograms in the data. If we actually
#kill these, we'll need to make template interferrograms that
#are at the location of the chopped interferrograms. (ie sidebands of chopper)


#create the template interferrogram:
  res = 1.0*30/(2*length)  #GHz spacing of points in the interferrogram, let's try to match resolution of the data
  intf = create_interf([],[],band=band, res=res)

##  if template:
##    restore, template
##    intf.x = tempx
##    intf.intf = tempy


  #now need to get the right number of pts to match the data:
  pt_space = optv/(1.0*rate)
  netscale = (intf['x'][1] - intf['x'][0])/pt_space
  xnew = pl.arange(pl.ceil(netscale*len(intf['x'])))*pt_space + min(intf['x'])
  intf2 = pl.interp(xnew, intf['x'],intf['intf'])
  regcut = .7    #tweak this around (?)

  q = pl.where(abs(xnew) < regcut)   #how far out in x to take the interferrogram, make this a function of frequency
  intf2 = intf2[q]
  intf2 = intf2.real


  #Okay, now we need to find the peaks:

  conv = find_intf(adc_cts2, intf2)

  #check that there is some S/N in the convolution function:
  convs_n = max(conv['conv'])/((conv['conv']).std())
  if convs_n < 4:
    outerr = {'scan_length':-1}
    print "too low S/N on conv " + str(convs_n)
    pl.figure()
    pl.subplot(2,1,1)
    pl.plot(pl.arange(len(conv['conv'])),conv['conv'])
    pl.title('Convolution Function for bolo ' + str(bolo))
    pl.subplot(2,1,2)
    pl.plot(pl.arange(len(adc_cts2)), adc_cts2)
    pl.title('Time Stream')
    return outerr


  print 'Conv S/N %s'%(convs_n)

  out = find_fts_peaks_max(conv, optv, rate, thres, length, adc_cts2,onepk=onepk, nscans=nscans)
  cpks = out['cpks']
  tpks = out['tpks']

  cpks.sort()

#-----------------------------------
  #Plot some stuff:

  pl.figure()
 
  pl.subplot(2,1,1)  
  xt =  (1.0*optv/rate)*pl.arange(len(adc_cts2))
  pl.plot(xt,adc_cts2,linewidth=0.5)
  pl.xlim(min((1.0*optv/rate)*conv['index']), max((1.0*optv/rate)*conv['index']))
  pl.ylim(min(adc_cts2), max(adc_cts2))
  pl.title('Time Stream')
  pl.xlabel('Optical Delay (cm)')
  if len(tpks) > 1:
    pl.plot(xt[tpks], adc_cts2[tpks], 'D',markerfacecolor='none',markeredgecolor='r')
  if len(tpks) == 1:
    tpks = pl.concatenate((tpks, tpks))
    pl.plot(xt[tpks], adc_cts2[tpks], 'D',markerfacecolor='none',markeredgecolor='r')
  pl.subplot(2,1,2)
  pl.plot((1.0*optv/rate)*conv['index'], conv['conv'],linewidth=0.2)
  pl.xlim(min((1.0*optv/rate)*conv['index']), max((1.0*optv/rate)*conv['index']))
  pl.title('Conv. function')
  pl.xlabel('Optical Delay (cm)')
  pl.plot(xt[cpks], conv['conv'][cpks], 'D',markerfacecolor='none',markeredgecolor='r')

  if len(include) != 0:
    cpks = cpks[include -1.0]


  if len(exclude) != 0:
    exclude = exclude - 1.0
    tind = pl.arange(len(cpks))
    for i in range(len(exclude)):
      n = where(tind != exclude[i])
      tind = tind[n]
      cpks = cpks[n]


#Now take the FFT and fill in the outputs
  fft_data = fts_fft(adc_cts2, optv, rate, cpks, length, bolo=bolo, band=pband, hann=hann, plt=plot, absv=absv, phaseb=phaseb,chop=chop1, crange=chrange, notquiet=notquiet)


  return fft_data




def time_fft(data2, samplerate=100., inverse=False,hann=False):
  '''
  IF N_PARAMS() EQ 0 then begin
     print, 'time_fft, data, samplerate=samplerate, inverse=inverse,hann=hann'
     return, -1
  ENDIF
  '''
  data=data2

  if hann:
    w1 = pl.hanning(len(data))
    data = data*w1

  #frequency axis:
  freqs = pl.arange(1+len(data)/2)/float(len(data)/2.0)*samplerate/2.                  # wut.
  if len(data) % 2 == 0 : freqs = pl.concatenate((freqs, -freqs[1:(len(freqs)-1)][::-1]))
  if len(data) % 2 != 0 : freqs = pl.concatenate((freqs, -freqs[1:len(freqs)][::-1]))

  response = pl.fft(data)
  if inverse : response = pl.ifft(data)

  out = {'freq': freqs, 'real': response.real, 'im': response.imag, 'abs': abs(response)}



  return out




def add_fts_scan(fts, ints=[], length=[]):

  '''
  IF N_PARAMS() EQ 0 then begin
    print, 'out = add_fts(fft_structure, int=int)'
    print, 'where /int adds the interferrograms, not the spectra'
    return, 0
  ENDIF
  '''
  if len(length) == 0 : length = 0
  if length == 0 : t = pl.where(fts['scan_length'] == max(fts['scan_length']))[0]
  if length != 0 : t = pl.where(fts['fft_length'] == length)[0]
  nscans = len(t)

  if len(ints) == 0:   #add the spectra, not the interferrograms. 
    freq = fts['freq'][:,t[0]]
    real = fts['real_FFT'][:,t[0]]
    im = fts['im_FFT'][:,t[0]]
    stdev = pl.zeros(len(freq))

    if nscans > 1:
      for i in range(1,len(t)):
        real = real + fts['real_FFT'][:,t[i]]
        im = im + fts['im_FFT'][:,t[i]]

      for j in range(len(freq)):
        stdev[j] = (fts['real_FFT'][j,t]).std()/pl.sqrt(nscans-1)
        #std dev of mean: std_Dev/sqrt(N-1)


    out = {'freq': freq, 'real': real/nscans, 'im': im/nscans, 'stdev': stdev}

  else:
    intf = fts['intf'][t[0]:]
    stdev = pl.zeros(len(intf))
    if nscans > 1:
      for i in range(1,len(t)):
        intf = intf + fts['intf'][:,t[i]]
      for j in range(len(intf)):
        stdev[j] = (fts['intf'][j,t]).std()/pl.sqrt(nscans-1)

    out = {'xint': reform(fts['xint'][t[0]:]), 'intf': reform(intf)/nscans, 'stdev': stdev}


  return out
