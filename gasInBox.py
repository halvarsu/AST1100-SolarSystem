import numpy as np
import random as ran
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt


N = 100000
k = 1.38864852e-23
m = 2*1.6737236e-27
T = 10000 # K




mu = 0
box_size = 1e-6
hole_size = box_size / 2
hole = (box_size/2-hole_size/2, box_size/2+hole_size/2)
dim = 3

seed = 27
np.random.seed(seed)

sigma = np.sqrt(k*T/m)
pos = np.random.uniform(0, box_size, (int(N), dim)) 
vel = np.random.normal(mu,sigma,(int(N), dim))
abs_vel = np.linalg.norm(vel, axis = 1)
x_vel = vel.T[0]
y_vel = vel.T[1]


start = 0
stop = 1e-9
n = 1000
dt = (stop-start)/n


#for k in xrange(n):

momentum_neg_z = 0
outside_count = 0
lower_wall_hit_count = 0
gained_momentum = 0
hit_hole = 0

for t in xrange(n):
    percent = t/float(n)*100
    status = "%2.2f%%"  % percent
    status = status + '\b'*(len(status)+ 1)
    print status,

    mask1 = pos > box_size # bool array, true if pos larger than box_size
    mask2 = pos < 0 

    # Summing momentum of the wall hitting
    momentum_neg_z += np.sum(np.abs(vel.T[2][mask2.T[2]])) * 2*m
    lower_wall_hit_count += np.sum(mask2.T[2])
            
    # Ma(s)king a hole in the box
    # x- and y-dimensions
    x_mask = np.logical_and(pos.T[0]>hole[0] , pos.T[0]<hole[1])
    y_mask = np.logical_and(pos.T[1]>hole[0] , pos.T[1]<hole[1])
    z_mask = pos.T[2] < 0

    # combining masks
    xy_mask = np.logical_and(x_mask, y_mask)
    xyz_mask = np.logical_and(xy_mask, z_mask)

    # summing
    hole_hit = np.sum(xyz_mask)
    outside_count += hole_hit
    gained_momentum += np.sum(np.abs(vel.T[2][xyz_mask])) * m

    vel[mask1] *= -1     # changes the indexes of vel where mask is true
    vel[mask2] *= -1     # see above
    accurate = False
    if accurate:
        pos = np.random.uniform(0,box_size,(N, 3)) 
        vel = np.random.normal(mu,sigma,(N, 3)) 
    #pos[xyz_mask] = np.random.uniform(0,box_size,(hole_hit, 3)) 
    #vel[xyz_mask] = np.random.normal(mu,sigma,(hole_hit, 3)) 
    pos += vel*dt


force_on_lower_wall = momentum_neg_z / (stop-start)
force_up = gained_momentum / (stop-start)

force_wanted = 914000.
boxes_needed = force_wanted/force_up
area_needed = hole_size**2 * boxes_needed

print ""
print "--------------------------------------"
print "    Number of particles: ", N
print "    Temperature: %d K"% T
print "    k:  ", k
print "--------------------------------------"
print "    Momentum in negative z-direction: ", momentum_neg_z
print "    Force on lower wall : " , force_on_lower_wall
print "    Hitting lower wall: ", lower_wall_hit_count
print "    Escaping: " , outside_count
print "    Upward force per box :    %g N" %force_up
print "    Boxes needed for %3dN:    %g"%(force_wanted , boxes_needed)
print "    area needed for  %3dN:    %g"%(force_wanted , area_needed)
print "    sides of area:               %g"% np.sqrt(area_needed)
print "    l_f/(4*u_f) = ", force_on_lower_wall/(force_up*4)


abs_vel2 = np.linalg.norm(vel, axis = 1)
kinetic_comp = 0.5* m * np.sum(abs_vel**2)/N
kinetic_comp2 = 0.5* m * np.sum(abs_vel2**2)/N
kinetic_calc = (3./2)*k*T

print "--------------------------------------"
print "    Computed kinetic energy:    ", kinetic_comp 
print "    Computed after sim:         ", kinetic_comp2
print "    Calculated kinetic energy:  ", kinetic_calc 
print "    Relative diff:              ",\
                    abs(kinetic_comp - kinetic_calc)/kinetic_comp
print "--------------------------------------"

plot = False
if plot:
    fig= plt.figure()
    num_bins = 50
    n_1,bins_1, patches = plt.hist(abs_vel, num_bins, alpha=0.5)
    n_2, bins_2, patches = plt.hist(abs_vel2, num_bins, alpha=0.5)
    print n_1, n_2
    print bins_1, bins_2
    plt.legend(['before', 'after'])
    s = 0
    for x in n_1:
        s+= x
    print s
    s=0
    for x in n_2:
        s+= x
    print s
    plt.show()
