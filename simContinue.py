"""
File: simContinue.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Goes trough different values made by MySateliteSim in folder
data and continues them
"""

from MySateliteSim import MySateliteSim
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



system = MySateliteSim(87464)
au = system.au
start_launch_v = 9585
for x in range(11):
    launch_v = start_launch_v + x
    p = np.load('positionsHomePlanet.npy')
    t = np.load('times.npy')
    st = np.load('data/times%d.npy' %launch_v)
    s = np.load('data/pos%d.npy' % launch_v)

    posFunc = interp1d(t[:-1], p)
    planetPos = posFunc(st)

    satFunc = interp1d(st, s.T)
    def vel(t,dt=0.0000001):
        return (satFunc(t+dt)-satFunc(t-dt))/(2*dt)

    planet4 = planetPos[:,4]
    relPos = s.T-planet4
    print relPos
    relDist = np.linalg.norm(relPos, axis =0)
    closest_dist = np.amin(relDist[:-10])
    print closest_dist
    closest_i = np.argmin(relDist[:-10])
    print closest_i
    relposCl = relPos[:,closest_i]
    hCl = closest_dist - system.radius[4]*1000./au
    tCl = st[closest_i]
    pCl = s[closest_i]
    aCl = np.arctan(relposCl[1]/relposCl[0])#angle
    vCl = vel(tCl)

    if hCl < 0:
        print "Wont bother with values inside planet"
        continue

    print "time of closest approach: ", tCl
    print "satpos at closest approach: ", pCl
    print "distance from surface at closest approach: ", hCl
    print "relpos at closest approach", relposCl
    print "angle at closest approach", aCl
    
    boost = 0
    angle = 0
    sim_length = 6
    
    dt_far = 1/(50000.)
    print "---------------------------------------------"
    print "continuing the launch ", launch_v
    system.boost_init(pos0, vel0, t0, init_boost = boost, 
            boost_angle = angle)
    psim, vsim, asim, tsim = system.satelite_sim(
            t0 = tCl, tN = tCl+sim_length, start_planet=4, 
            target_planet =5, vel0 = vel0, pos0 = pos0, 
            first_boost = False, dt_far=dt_far, sec_per_dt_close = 5)
    filenames = ['%s%d' %(name,launch_v) for name in ('pos','vel','times')]
    system.saveData(fname = filenames, folder = 'postGAdata')
#system.plot()

