import numpy as np
from MySateliteSim import MySateliteSim
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


if __name__ == "__main__":
    system = MySateliteSim(87464)
    pos = np.load('data/orbit_pos.npy').T
    vel = np.load('data/orbit_vel.npy').T
    t = np.load('data/orbit_times.npy')
    tOrbit = np.load('data/tOrbit.npy')

    tStop = t[-2]
    print tStop
    posFunc = interp1d(t,pos)
    velFunc = interp1d(t,vel)
    #times = np.linspace(0, tOrbit, 10)
    times = np.linspace(0,tStop,10)
    pos = posFunc(0)
    vel = velFunc(0)
    r = np.linalg.norm(pos)
    v = np.linalg.norm(vel)
    print pos/r, vel/r, system.horizonVec(pos,vel)

    u_h = system.horizonVec(pos,vel)
    phi = np.arctan2(u_h[1], u_h[0])
    theta = np.arccos(u_h[2])
    print 'video ', 1.1, theta, phi, pos[0], pos[1], pos[2] 
    for i, t in enumerate(times):
        pos = posFunc(t)
        vel = velFunc(t)
        u_h = system.horizonVec(pos,vel)
        u_v = vel/np.linalg.norm(vel)
        r = np.linalg.norm(pos)
        #phi =  np.arctan2(-pos[1],-pos[0]) 
        #theta = np.arccos(-pos[2]/r)
        phi = np.arctan2(u_h[1], u_h[0])
        theta = np.arccos(u_h[2])
        print 'picture ', t, theta, phi, pos[0]/r, pos[1]/r, pos[2]/r

    print 'video ', tStop+0.1, theta, phi, pos[0]/r, pos[1]/r, pos[2]/r
