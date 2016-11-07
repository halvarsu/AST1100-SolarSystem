import matplotlib.pyplot as plt
import sys
from AST1100SolarSystemViewer import AST1100SolarSystemViewer
import argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years',default=200,
                    type=int,help='years to run simulation')
parser.add_argument('-r', '--res',default=4000,
                    type=int,help='timesteps per year')
parser.add_argument('-c', '--checkplanetpositions',default=False,
                    action='store_true', help=('test planet positions vs '\
                            '"exact" positions from AST1100SolarSystem'))
parser.add_argument('-s', '--save',default=False,action='store_true',
                    help='save times as times.npy')
parser.add_argument('-l', '--load',default=False,action='store_true',
                    help='load all data from data/temp*data.npy')
parser.add_argument('-i', '--integrator',default='euler_cromer',
                    choices=['leapfrog','euler_cromer'])
args = parser.parse_args()

years = args.years 
res   = args.res

from MySolarSystem import MySolarSystem
system = MySolarSystem(87464)

if args.load:
    pos,vel,accel = np.load('data/temp_posvelaccdata.npy')
    times = np.load('data/temp_timesdata.npy')
else:
    data = system.simulate(years = years, res = res, 
                               integrator = args.integrator)
    print "saving data/temp_data.npy"
    np.save('data/temp_posvelaccdata', data[:-1])
    np.save('data/temp_timesdata', data[-1])
    pos, vel, accel, times = data


linsp = np.linspace(0,1,years*res+1)

r = np.sqrt(pos[:,0]**2 +pos[:,1]**2)
theta = np.linspace(0,2*np.pi,res*years+1)


fig, ax = plt.subplots()
n = len(pos)
ax.plot(pos[:,0,:], pos[:,1,:])
ax.scatter(pos[:,0,:][::n/years],
        pos[:,1,:][::n/years])
ax.scatter(0,0, c='y')

plt.axis('equal')
plot_max = max(system.a +1)
plt.axis([-plot_max,plot_max,-plot_max,plot_max])
ax = system.polar_orbits(ax =ax)

planetPos = pos.swapaxes(0,1).swapaxes(1,2) 
testP = planetPos[:,:, :-1]

fig = plt.figure()
ax = plt.subplot(111)
plt.plot(accel[:,0], times)
plt.show()
if args.save:
    np.save('times', times)
    np.save('accel_data', accel)
if args.checkplanetpositions:
    system.checkPlanetPositions(testP, years, res)
