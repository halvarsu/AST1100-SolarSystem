"""
File: rocket_eq.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: 
"""
import numpy as np
import sys

def load_data(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    data = [int(lines[0])] # first line, outside_count
    for line in lines[1:]:
        data.append(float(line))
    return data

try:
    delta_v = float(sys.argv[1])
    #time =  float(sys.argv[2])
except IndexError:
    print "WrongUsageError! usage: python %s delta_v " %sys.argv[0]
    sys.exit()

print "Loading data from simulated engine..."
data = load_data('data_dump.dat')
outside_count   = data[0]
m               = data[1]
force_up        = data[2]
fuel_per_box    = data[3]
dt              = data[4]
gained_momentum = data[5]

v_e = gained_momentum/(m*outside_count)
force_needed = gained_momentum / dt
momentum_per_box = outside_count * m / dt

print "Calculating..."
# To be given: 
dry_mass = 1190
    
    

fuel_mass = dry_mass*np.exp(delta_v/v_e) - dry_mass

print "Satelite mass  : %11d kg"    % dry_mass 
print "Delta v needed : %11.2f m/s" % delta_v
print "Fuel needed    : %11.2f kg"  % fuel_mass
#print "Boxes needed    : %11.2f s"   % boxes




