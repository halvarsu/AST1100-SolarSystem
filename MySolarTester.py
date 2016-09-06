import matplotlib.pyplot as plt
from MySolarSystem import MySolarSystem

seed = 87464
system = MySolarSystem(seed)
years = 200; res = 10000
pos, accel = system.simulate(years = years, res = res)

import numpy as np
linsp = np.linspace(0,1,years*res+1)

r = np.sqrt(pos[:,0]**2 +pos[:,1]**2)
#plt.plot(linsp, accel[:,0,4])
#plt.plot(linsp, pos[:,0,4])
#plt.plot(linsp, r[:,4])
#plt.legend(['a','x'])
fig, ax = system.planet_plotter(pos, 10)
ax.scatter(0,0, c='y')
plt.axis('equal')
#plot_max = max(system.a +1)
#plt.axis([-plot_max,plot_max,-plot_max,plot_max])
#ax = system.polar_orbits(ax =ax)
plt.show()
