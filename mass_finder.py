"""
File: rocket_sim3.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Information copyed from simulation by gasInBox.py used to
calculate longer rocket burns for all your burning desires
"""

def mass_finder(boost, test = False, mass = 1100):
    import numpy as np
    import matplotlib.pyplot as plt
    from MySolarSystem import MySolarSystem, load_engine_data


    def ms_to_auyear(speed):
        au   = 149597871 #km
        year = 3600*24*365.25 #seconds
        return  speed*au*1000/year #m/s to au/year

    data = load_engine_data('gasInBox/data_dump.dat')
    outside_count   = data[0] # particles leaving per dt
    m               = data[1] # mass per particle
    force_up        = data[2] # force per box. = gained_momentum/dt
    fuel_per_box    = data[3] # outside count * m
    dt              = data[4] # timestep for simulation
    gained_momentum = data[5] # momentum per box per dt


    dpdt = gained_momentum/dt
    num_p = outside_count/dt

    # boost_time should be less than 20 min
    time = 60

    force_tot = boost * mass / time
    num_boxes = force_tot / dpdt
    fuel =  num_p * num_boxes * m  *time# /dt #fuel_per_box * num_boxes * time #

    size_per_box =  1e-6
    area = size_per_box **2 * num_boxes
    np.save('data_dump2', (5000, dpdt, num_boxes, num_p))
    

    if test:
        print "\nA burn time of %d seconds gives the following parameters:\n" %time
        print num_boxes, boost, num_p, dpdt
        print 'area of boxes: ', area ,'m (area of "nozzle" is 1/4 of this)'
        print "force needed: ", force_tot , 'N'
        print "fuel needed: ", fuel , "kg\n"

        seed = 87464
        system = MySolarSystem(seed)
        system.massNeededCheck(num_boxes, boost, dpdt,
            num_p, fuel)
        boost = ms_to_auyear(boost)
    return fuel

if __name__ == "__main__":

    import sys
    try:
        boost= int(sys.argv[1])
    except IndexError:
        print "Bad usage, {0} takes one argument.\nUsage: python {0} boost_wanted".format(sys.argv[0])
        sys.exit()
    mass_finder(boost, test = True)
