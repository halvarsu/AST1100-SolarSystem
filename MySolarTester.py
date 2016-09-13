import matplotlib.pyplot as plt
import sys
from MySolarSystem import MySolarSystem
from AST1100SolarSystemViewer import AST1100SolarSystemViewer

seed = 87464
system = MySolarSystem(seed)
years = 200; res = 4000
#print "CALCULATING LEAPFROG"
try:
    arg1 = sys.argv[1]  
except IndexError:
    arg1 = 'normal'
print "--- starting %s simulation ---" %arg1
if arg1 == 'nbody':
    pos, accel, times = system.simulate_nbody(years = years, res = res)
elif arg1 == 'normal':
    pos, accel, times = system.simulate(years = years, res = res, 
                        integrator = 'euler_cromer')
else:
    pos, accel, times = system.simulate(years = years, res = res)

import numpy as np
linsp = np.linspace(0,1,years*res+1)

r = np.sqrt(pos[:,0]**2 +pos[:,1]**2)
theta = np.linspace(0,2*np.pi,res*years+1)

fig, ax = system.planet_plotter(pos, 10)
ax.scatter(0,0, c='y')
plt.axis('equal')
plot_max = max(system.a +1)
plt.axis([-plot_max,plot_max,-plot_max,plot_max])
ax = system.polar_orbits(ax =ax)

print pos.shape
planetPos = pos.swapaxes(0,1).swapaxes(1,2) 
print planetPos.shape
testP = planetPos[:,:, :-1]
print testP.shape

system.checkPlanetPositions(testP, years, res)
#system = AST1100SolarSystemViewer(seed)
#system.orbitXml(planetPos, times)

fig = plt.figure()
ax = plt.subplot(111)
np.save('times', times)
np.save('accel_data', accel)
plt.plot(accel[:,0], times)
plt.show()
