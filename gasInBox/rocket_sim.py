"""
File: rocket_sim.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Information copyed from simulation by gasInBox.py used to
calculate longer rocket burns for all your burning desires
"""

import numpy as np
import matplotlib.pyplot as plt
from rocket_test import test_rocket

from load_data import load_data


data = load_data('data_dump.dat')
outside_count   = data[0]
m               = data[1]
force_up        = data[2]
fuel_per_box    = data[3]
dt              = data[4]


time_step = 0.1
force_wanted = 500
start_fuel = 688
boxes_needed = force_wanted/force_up
momentum_per_box = outside_count * m *time_step / dt

total_momentum = momentum_per_box * boxes_needed*time_step /dt
total_fuel_use = fuel_per_box * boxes_needed*time_step /dt

sat_mass = 1190
fuel = start_fuel
mass = sat_mass + fuel

a = force_wanted / float(mass)
dv = 0
time = 0
print 'm ', m            
print 'momentum per box =', force_up     
print 'fuel_use_total   =', total_fuel_use 
print 'boxes_needed     =', boxes_needed
print dt           

while fuel > 0:
    status = "fuel left: %2.2f%%   dv total: %g"  % (fuel, dv)
    status = status + '\b'*(len(status)+ 1)
    print status,

    
    a = force_wanted / float(mass)
    dv = dv + a *time_step
    fuel -= total_fuel_use
    time += time_step

print ""
print "time spent: %d s = %d min" %(time, time/60)
print "delta v   : ", dv
print ""
    
test_rocket(boxes_needed, dv, momentum_per_box, outside_count/dt, start_fuel)
