import numpy as np
import sys
import matplotlib.pyplot as plt
from MySateliteSim import MySateliteSim
from rocket_eq import load_data

system = MySateliteSim(87464)

t0 = 14.319907707
system.launch_init(t0 = t0, init_boost = 7338.0869, 
                   init_angle = 0.6*np.pi, start_dist = 0,
                   ground_launch = False)
write = False
if write:
    with open('launch_text.txt','a') as outfile:
        planet_vel0 = system.velFunc(t0)[:,0]
        outfile.write('delete this\n')
        engine_specs = np.load('data_dump2.npy')
        print engine_specs
        x,y = system.pos0
        vx, vy = system.vel0 - planet_vel0
        fuel, dpdt, num_boxes, num_p = engine_specs
        data =  [t0, x, y, vx, vy, 5000, dpdt, num_boxes, num_p]
        data_str = [str(d) for d in data]
        launch_line = 'launch ' + ' '.join(data_str)
        print launch_line
        outfile.write(launch_line)


pos1, vel1, acc1, times1 = system.satelite_sim(
        tN = 14.32-t0, speed_factor = [0.01, 1, 1], target_planet = 0,
        break_close = False)

planet_vel = system.velFunc(times1[-1])[:,0].T
orbit_vel  = vel1[-1] - planet_vel
orbit_dir  = orbit_vel/np.linalg.norm(orbit_vel)
orbit_boosted = system.boost(orbit_vel, 3337.9885, boost_dir = orbit_dir)
new_vel = orbit_boosted + planet_vel
print 'prev vel : ', vel1[-1]
print 'new vel  : ', new_vel
print 'delta vel: ', new_vel - vel1[-1]
print 'last pos : ', pos1[-1]
print 'last time : ', times1[-1]


print orbit_boosted - orbit_vel

raw_input()
# For circulization: 264.91

system.boost_init(pos0 = pos1[-1], vel0 = new_vel, t0 = times1[-1], 
                  init_boost = 0)

pos2, vel2, acc2, times2 = system.satelite_sim(
        tN = 3.0004, speed_factor = [0.5, 1, 1], target_planet = 4)

pos = np.concatenate((pos1[:-1],pos2[:]))
vel = np.concatenate((vel1[:-1],vel2[:]))
times = np.concatenate((times1[:-1],times2[:]))

#planet4 = system.posFunc(times)[:,4]
#relpos4 = planet4 - pos.T
#km_relpos4 = 1./1000 * system.au * relpos4
#print km_relpos.T[-1], km_relpos.T[0]
#print "Distance to planet surf(km): ", \
        #system.radius[0] - np.linalg.norm(km_relpos.T[-1])

system.plot(planet=4)

planet4 = system.posFunc(times)[:,4]
relpos4 = planet4 - pos.T

print 'last time: ', times[-1]
rel_r4 = np.linalg.norm(relpos4,axis = 0)
closest_i = np.argmin(rel_r4)
print 'time of closest', times[closest_i]
if raw_input('save?(y/N)') == 'y':
    b = 'fromGround'
    fnames = [w+b for w in ['pos','vel','times']]
    system.saveData(fname = fnames)
