"""
File: rocket_sim.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Information copyed from simulation by gasInBox.py used to
calculate longer rocket burns for all your burning desires
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from rocket_test import test_rocket

from load_data import load_data


data = load_data('data_dump.dat')
outside_count   = data[0] # particles leaving per dt
m               = data[1] # mass of particles
force_up        = data[2] # per box. = gained_momentum/dt
#fuel_per_box    = data[3] # outside count * m
dt              = data[4] # timestep for simulation
gained_momentum = data[5] # momentum from all boxes

try:

    boost_wanted = int(sys.argv[1])
except IndexError:
    print "Bad usage, {0} takes one argument.\nUsage: python {0} boost_wanted".format(sys.argv[0])
    sys.exit()


momentum_per_box = gained_momentum/ force_up/dt
mass = 1190 #kg

time = 19*60
force_needed = mass*boost_wanted/time
boxes_needed = force_needed/force_up
fuel_needed = outside_count * m * time/dt  *boxes_needed

print boxes_needed
print fuel_needed
#a = force_up / mass

test_rocket(boxes_needed, boost_wanted, momentum_per_box, outside_count/dt,
        fuel_needed)
