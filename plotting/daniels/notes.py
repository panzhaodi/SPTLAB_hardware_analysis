tuner = ppf.DeviceTuner(devicesDirectory = hardware_map,
                        verbose = True,
                        saveData = True,
                        outputDirectory = data_directory,
                        is16Channel = True, initializeBoards = False)
ppf.initializePPTunerForUse(tuner)

resies = ppf.overbiasAndNull(tuner, squids = sSquids, amplitude = BP.default_overbias_amp, channels = 'used',
                                 carrier_gain = BP.default_carrier_gain, nuller_gain = BP.default_nuller_gain, 
                                 excludedChannels = BP.unbiasable_channels, specificArgs = BP.overbias_specargs, 
                                 verbose =True,fast=True, leave_dan_on=True,units='NORMALIZED_16ch')#useStoredValues=True



    fDic = getFrequencyMapPPChanForm(tuner)

def getFrequencyMapPPChanForm(tuner):
    freqMap = {}
    for squidId in tuner.squidIds:
        squid = tuner.deviceManager.findDevice(SQUID_TYPE, squidId)
        freqs = squid.getCombFrequencies()
        nchans = dfmux_conf['NumDMFSChans']
        if type(freqs)!=type({}):
            print 'PyPolInterface.getFrequencyMapPPChanForm could not find any bias frequencies in the hardware map for squid %s, using zeroes.' % squidId
            for c in range(1, 1+nchans):
                channel_name = 'Ch%i' % c
                freqMap['%s%s'%(squidId,channel_name)] = {'frequency':0}
        else:
            for c in range(1, 1+nchans):
                channel_name = 'Ch%i' % c
                if channel_name in freqs.keys(): freq = freqs[channel_name]
                else: freq = 0 # the channels that don't have associated frequencies get zeroes
                freqMap['%s%s'%(squidId,channel_name)] = {'frequency':int(float(freq))}
    #print freqMap
    return freqMap


ppCallTF(tuner, 'set_frequency',  squids=squids, channels=channels, specificArgs = fDic, excludedChannels = excludedChannels, 
             frequency = 0, units = 'Hz', target = 'carrier')

def ppCallTF(tuner, funcName, boards = [], squids=[], channels=None, sqChs = [],  specificArgs = {}, excludedChannels = [], **kargs):
    '''
    Runs tuber functions in parallel

    tuner: ppTuner
    funcName: the name of the tuber function
    specificArgs:  see DFMLDoc in root directory

    Include the arguments to the tuber function as keyword arguments to this function  (this is what the **kargs corresponds to)

    for all other arguments see ppGetAddressList

    returns dic of the form {ppAddress: results of tuber func}
    '''


    callTarget = ppGetCallTarget(tuner, boards, squids, channels, sqChs = sqChs, excludedChannels = excludedChannels)
    return callTF( tuner.DFML_bmm, callTarget, funcName, specificArgs = specificArgs, addressConverter = tuner.DFML_ac, **kargs)

def ppGetCallTarget(tuner, boards = [], sqLst = [], chLst = None, sqChs = [], excludedChannels = []):    
    '''
    tuner: ppTuner
    for the other arguments see the ppGetAddressList documentation

    returns: a callTarget for devices specified
    '''
    callTarget = addressListToCallTarget(ppGetAddressList(tuner, boards = boards, sqLst = sqLst, chLst = chLst, sqChs=sqChs, excludedChannels = excludedChannels))
    for k in tuner.DFML_DeadBoards:
        if k in callTarget:
            callTarget.pop(k)
    return callTarget

#if chLst is all it will fill in all the channels
def ppGetAddressList(tuner, boards = [], sqLst = [], chLst = 'all', sqChs = [], excludedChannels = []):    
    '''
    There are lots of options for returning the devices we want to work with.  Sorry if this isn't clear.


    Returns: a list of the DFML style addresses for the devices specified.  DFML style address is the board#_mux#_chan# format so '90_1_3' would be valid

    boards: the ids of the boards we want to use if we only care about interacting with the boards. If boards != [] all the following arguements are ignored.
    sqLst:  A list of the squids we are working on.  If we wish to use all the squids in the tuner, set sqLst to the string "all".
    chLst:  A list of the channels we are working on.  If chLst == "all", we use channels 1-16.  If chLst == 'used' we only use the channels that have stored 
            frequencies in the hardware map.  If chLst = None (where none is the python object None, not the string) we only a list of the addresses for the squids.


    excludedChannels:  if we wish to exclude specific channels on specific squids we can specify it through this.  excludedChannels is a list of pypol style addresses
                       that specify the squids/channels we do not included.  

    sqChs: If we wish to add specific channels on specific squids we specify this with sqChs.   This is a list of the pypol style addresses (e.g. 'Sq3Ch2')
    
    '''
        

    if len(sqLst) == 0 and len(sqChs) == 0:
        if boards == 'all':
            return tuner.DFML_bmm.keys()
        if len(boards) == 0:
            raise RuntimeError("ppCallTF should be called on some board or some squid")
        else:
            return boards

    if sqLst == 'all':
        sqLst = tuner.squidIds

    if type(sqLst)!=type([]):
      raise Exception('WARNING: PyPolInterface.ppGetAddressList was provided with non-list value for the sqLst argument: %s' % str(sqLst))
    
    addressLst = []

    for sqch in sqChs:
        if sqch not in  excludedChannels:
            addressLst.append(tuner.DFML_ac.getAddress(sqch[:sqch.rfind('Ch')])+'_'+sqch[sqch.rfind('Ch')+2:])


    for sq in sqLst:
        if chLst == None:
            addressLst.append(tuner.DFML_ac.getAddress(sq))
            continue

        if chLst == 'all':
            localChLst = ['Ch%d'%c for c in range(1, dfmux_conf['NumDMFSChans']+1)]
        elif len(chLst) == 0 or chLst == 'used':
            localChLst = [str(k) for k in tuner.get_squid(sq).getCombFrequencies()]
        elif chLst =='first':
            localChLst = [str(k) for k in tuner.get_squid(sq).getCombFrequencies()][0:1]
        elif chLst =='fan':
            localChLst = [str(k) for k in tuner.get_squid(sq).getCombFrequencies()][0:2]
        elif chLst =='second':
            localChLst = [str(k) for k in tuner.get_squid(sq).getCombFrequencies()][1:2]
        else:
            localChLst = chLst
        for ch in localChLst:
            if sq+ch not in excludedChannels:
                addressLst.append(tuner.DFML_ac.getAddress(sq,ch))
    return addressLst





