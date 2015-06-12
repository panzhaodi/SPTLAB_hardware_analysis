import pyEagle.Eagle as Eagle
import pyEagle.Utils as Utils
import math

SIGNAL_LAYER = 1

# And Pitches
wirebond_pitch = 0.1 
zif_pitch = 0.3

cable = Eagle.Board(outfile='flexCable.scr')

wirebond = Eagle.FootprintFactory("WIREBOND-TEST", "WIREBOND-NEW")
zif1 = Eagle.FootprintFactory("ZIF-NEW", "ZIF-NEW1")
zif2 = Eagle.FootprintFactory("ZIF-NEW", "ZIF-NEW2")

def dsin(deg):
    """Sin in degrees."""
    return math.sin(math.radians(deg))

def dcos(deg):
    """Cos in degrees."""
    return math.cos(math.radians(deg))

wafer = wirebond("Wafer","R0",(0,-1.2))
cable.add(wafer)

v_spacing = 22  # Placement of the ZIF connectors.
h_spacing = 15

zif_base=85

bolos_per_zif=33  # 


connector_angles=[-60,-45,-34,-0,0,25,45,70]
connector_start_mods=[-1.2,2.8,6,0,0,7.8,5.7,2.8]
connector_lengths=[30,25,20,5,5,20,20,30]


zif_xstarts=[39.7,45.8,47.4,33.0,33.0,47.5,40.5,32.9]
zif_ystarts=[-24.0,-8.6,5.4,23.4,30.6,45.0,57.2,76.3]

for board in xrange(4): # Branch off for each LC board.

  for connector in [0,1]: # Then branch off for each zif connector.

    if board==0:
      x = zif_base + h_spacing*(connector % 2)
    if board==1:
      x = zif_base + h_spacing*(connector % 2)
    if board==2:
      x = zif_base + h_spacing*(((1+connector) % 2)+1)
    if board==3:
      x = zif_base - h_spacing*(connector % 2)

    y = -63 + v_spacing*(2*board+connector)
    if connector == 0 :       
      zif_pad = zif1("ZIF_CON%d" % (2*board+connector), "R0", (x, y))
    if connector == 1 :
      zif_pad = zif2("ZIF_CON%d" % (2*board+connector), "R0", (x, y))
    cable.add(zif_pad)

    connector_angle=math.pi*connector_angles[2*board+connector]/180


    zifstart_x = zif_xstarts[2*board+connector]
    zifend_x = x-20
    zifstart_y = zif_ystarts[2*board+connector]
    zifend_y = y+zif_pitch*44

    zif_angle = math.atan2(zifend_y-zifstart_y,zifend_x-zifstart_x)

    if 2*board+connector != 3:
      if 2*board+connector !=4:
        zif_angle=connector_angle
    if 2*board+connector != 4:
      if 2*board+connector !=3:
        zif_angle=connector_angle


    for pad in xrange(2*bolos_per_zif):
      pad_num = (2*board+connector)*2*bolos_per_zif+pad

      trace1 = Eagle.Signal(0.04, SIGNAL_LAYER, 2) # Defining the thinner trace
      trace1.add(-.08, (pad_num+12)*wirebond_pitch-1.2)  # Adding the pads to signal
      trace1.r_theta(0.55,0) # Thin straight segment

      trace2 = Eagle.Signal(.045, SIGNAL_LAYER, 2) # Defining the thicker trace
      trace2.add(0.27, (pad_num+12)*wirebond_pitch-1.2)  # Starting this signal where the thin one left off.
      trace2.r_theta(1.3,0)

      y_shift = 6*wirebond_pitch
      trace2.add(3.5,trace2.last.y+y_shift)

      
      trace3 = Eagle.Signal(.05, SIGNAL_LAYER, 2)
      trace3.add(trace2.last.x,trace2.last.y)
      trace3.add(13,trace3.last.y)

      trace= Eagle.Signal(0.045,SIGNAL_LAYER,2)
      trace.add(trace3.last.x,trace3.last.y)
####### Traces are 45 um wide, 100 um apart. Have extra 5 um spacing that I can compress ######
####### I want the outer boards to compress to the outside edges, and the inner boards to compress to a line ~2/3 their width away from the center ###### 
      if board == 0:
        y_shift=-1*(connector*(2*bolos_per_zif-1)+pad)
      elif board == 3:
        y_shift=4*bolos_per_zif-2-(connector*(2*bolos_per_zif-1)+pad)
      elif board==1:
        y_shift=math.copysign(math.ceil(abs((4*bolos_per_zif-1)/3-(connector*(2*bolos_per_zif-1)+pad))),(4*bolos_per_zif-2)/3-(connector*(2*bolos_per_zif-1)+pad))
      elif board==2:
        y_shift=math.copysign(math.floor(abs((4*bolos_per_zif-1)/1.5-(connector*(2*bolos_per_zif-1)+pad))),(4*bolos_per_zif-2)/1.5-(connector*(2*bolos_per_zif-1)+pad))

      squeeze_angle=math.atan2(0.005*y_shift,3)  # I want this compression to occur over ~ 3mm
## Set up squeezes ##
      if board==0:
        counter=connector*(2*bolos_per_zif-1)+pad # starts at 0 and goes to max
      elif board==3:
        counter=(4*bolos_per_zif-2)-(connector*(2*bolos_per_zif-1)+pad) # starts at max and goes to 0
      else:
        counter=math.ceil(abs((4*bolos_per_zif-2)/2.-(connector*(2*bolos_per_zif-1)+pad))) # starts and stops at max, 0 in middle

      trace.r_theta(counter*abs(0.10/math.tan((math.pi-squeeze_angle)/2)), 0) # set up squeezes
      
## Do the squeezes ##
      if squeeze_angle != 0:
        trace.r_theta(.005*y_shift/math.sin(squeeze_angle),squeeze_angle)
## Now move on ##
      trace.add(20, trace.last.y)

      trace.r_theta(2+connector_start_mods[2*board+connector]+pad*0.10/math.tan((math.pi+connector_angle)/2.),0) #set up connector split

      trace.r_theta(5+pad*0.10/math.tan((math.pi+connector_angle)/2.),connector_angle) # connector split
 
#      trace.r_theta(connector_lengths[2*board+connector]+pad*0.10/math.tan((math.pi+connector_angle)/2.),connector_angle) # connector split

      if pad==33:
        print (trace.last.x,trace.last.y)


##### I want the zif splits to start off 45 um wide with 50 um spacing, but quickly increase to 50 um wide, 50 um spacing #####
      trace.r_theta(5-pad*0.10/math.tan((math.pi+(connector_angle-zif_angle))/2.),connector_angle) #set up zif split
      trace.r_theta(1,zif_angle)

      counter = math.copysign(math.floor(abs(pad-(2*bolos_per_zif-1)/2.)),pad-(2*bolos_per_zif-1)/2.)
      perp=.005*counter+math.copysign(1,counter)*0.0025
      parallel=5
      trace.r_theta(math.sqrt(perp**2+parallel**2), zif_angle+math.copysign(1,perp)*math.pi/2.-math.atan(parallel/perp))

      trace4 = Eagle.Signal(.05, SIGNAL_LAYER, 2)
      trace4.add(trace.last.x,trace.last.y)

      y_distance=trace4.last.y-(y+zif_pitch*44-(33.5-pad)*wirebond_pitch)
      x_distance=y_distance/math.tan(zif_angle)

      trace4.r_theta(math.sqrt(y_distance**2+x_distance**2),zif_angle) #zif split
      trace4.r_theta(.5,0)    

      trace4.add(x-15,trace4.last.y)

      cable.add(trace1)
      cable.add(trace2)
      cable.add(trace3)
      cable.add(trace)
      cable.add(trace4)

      if connector == 0:
        trace4.add(x-5-1.25*abs(trace4.last.y-y-14.7),trace4.last.y)
        trace4.add(x-3.5,y+pad*zif_pitch)
      else:
        trace4.add(x-5-1.25*abs(trace4.last.y-y-11.1),trace4.last.y)
        trace4.add(x-3.5,y+(90-2*bolos_per_zif+pad)*zif_pitch)
      trace4.r_theta(1.75*(1+pad % 2),0)


cable.draw()




