"""
File: simContinueFineTuning.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Used for finding the leg from gravity assist planet to target
planet, with different values for boost
"""

from MySateliteSim import MySateliteSim
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--part', choices=range(1,7),default=6,type=int)
args = parser.parse_args()

system = MySateliteSim(87464)
au = system.au


posFunc = system.getPositionsFunction()

if args.part == 1:
    archive_txt = 'last_launch.txt'
    Position = (3.60039137872,-0.551200231282)
    Velocity = (1.68787359777, 3.91496135905)
    t0 = 15.224431483
    tN = 1.597
    start_planet = -1
    target = 4
    sf = [1,2,4]
    boost_v = 30; boost_a = np.pi*0.9
    archive_txt = 'last_launch.txt'

    #b80 pi 0.00448145217402
    #b70 pi 0.00219013469487 
    #b60 pi 0.000634121339111 # Too much push out of system
    #b65 too much push

elif args.part == 2:
    archive_txt = 'last_launch2.txt'
    #Position = ( 2.01364712592 ,  4.30981258814  )
    #Velocity = ( -2.28924475191 ,  4.00828391069  )
    #t0 = 16.801431483

    t0 = 16.7928878928 
    Position = ( 2.04224112167 ,  4.29694107731  )
    Velocity = ( -2.80886338727 ,  -0.669233649925  )

    tN= 5.45
    start_planet = 4
    target = 5
    sf = [1,2,1]
    boost_v = 645.5; boost_a = 0
    #b580 a=0 0.27342445429
    #b600 0.191132077508
    #b645 0.00857565471599
    #b655 0.0297313401617
    #b649 0.00678046150003
    #b647 0.0106440945962
    # increasing accuracy
    #b647 0.00263887587696

elif args.part == 3:
    archive_txt = 'last_launch3.txt'
    Position = ( -2.19852605635 , -4.71012963539  )
    Velocity = ( 2.6946849101 , -2.28194618607  )
    t0 = 20.0660134278

    tN = 1.5
    start_planet = -1
    target = 5
    sf = [1,0.5,1]
    boost_v = 45; boost_a = np.pi/2
    #boost_v = 0; boost_a = 0
    #v40 a pi/2 0.000607432333933
    #v45 0.000121742333822               
elif args.part==4:
    t0 =  20.8568550951
    Position = ( 0.125152282389 , -5.87411880252  )
    Velocity = ( 3.03327034075 , -0.697120054535  )

    archive_txt = 'last_launch4.txt'
    tN = 1.4
    start_planet = -1
    target = 5
    sf = [1,0.5,1]
    boost_v = 0; boost_a =0

elif args.part == 5:
    archive_txt = 'last_launch5.txt'
    Position =  (1.68932738, -5.98807047)
    Velocity =  (3.25026227, -0.17369137)
    t0 = 21.3817592525

    tN = 0.015/16.05 #0.0099065#0.015/8.025#21.3818426335 - t0 
    start_planet = 5 
    target = 5
    sf = [1,1,1]
    boost_v = 2407.03
    boost_a = np.pi

elif args.part == 6:
    archive_txt = 'last_launch6.txt'
    t0 = 21.3817592526 
    Position = ( 1.68932736458 ,  -5.98807045114  )
    Velocity = ( 3.16245857718 ,  0.327102058371  )

    tN = 0.015/16.05 #0.0099065#0.015/8.025#21.3818426335 - t0 
    start_planet = 5 
    target = 5
    sf = [1,1,1]
    boost_v = 0
    boost_a = np.pi

pos0 = np.array(Position)
vel0 = np.array(Velocity)
if boost_v != 0:
    if start_planet == -1:
        new_vel = system.boost(vel0, boost_v, boost_a)
    else:
        velFunc = system.getVelocityFunction()
        planet_vel = velFunc(t0)[:,start_planet].T
        orbit_vel  = vel0 - planet_vel
        orbit_dir  = orbit_vel/np.linalg.norm(orbit_vel)
        orbit_boosted = system.boost(orbit_vel, boost_v, boost_angle = boost_a, boost_dir = orbit_dir)
        new_vel = orbit_boosted + planet_vel
        print 'pla vel: ', planet_vel
    print 'BOOST'
    print 'old vel: ', vel0
    print 'new_vel: ', new_vel
    delvel = new_vel - vel0
    print 'deltvel: ', delvel
    vel0 = new_vel
else:
    delvel = 0
print "---------------------------------------------"

if start_planet == target:
    eacp = False
else:
    eacp = True
system.boost_init(pos0 = pos0, vel0 = vel0, t0 = t0)
pos1, vel1, acc1, times1 = system.satelite_sim(
        tN = tN, speed_factor = sf, target_planet = target,
        break_close = eacp)

np.save('test_pos', pos1)
np.save('test_times', times1)
cl, cl_i = system.plot(planet = target, end_at_closest_approach=eacp,
        plot_length = -1)
print times1[cl_i], times1[-1]
with open(archive_txt,'w') as outfile:
    outfile.write('boost: ')
    outfile.write('\ndv: ' + str(delvel))
    outfile.write('\npos: '+ str(pos1[0]))
    outfile.write('\nvel: '+ str(vel1[0]))
    outfile.write('\nacc: '+ str(acc1[0]))
    outfile.write('\ntime: '+ str(times1[0]))
    outfile.write('\nAt closest approach: ' )
    outfile.write('\npos: '+ str(pos1[cl_i]))
    outfile.write('\nvel: '+ str(vel1[cl_i]))
    outfile.write('\nacc: '+ str(acc1[cl_i]))
    outfile.write('\ntime: '+ str(times1[cl_i]))
    outfile.write('\ndist: '+ str(cl))
print 'Wrote data to ', archive_txt
