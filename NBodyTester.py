import matplotlib.pyplot as plt
import numpy as np
import sys
from NBodySystem import NBodySystem

seed = 87464
system = NBodySystem(seed)
years = 100; res = 1000

pos, accel, times = system.simulate_nbody(years = years, res = res)

linsp = np.linspace(0,1,years*res+1)

#fig, ax = system.planet_plotter(pos, 10)
plt.axis('equal')
plot_max = max(system.a +1)
plt.axis([-plot_max,plot_max,-plot_max,plot_max])
#ax = system.polar_orbits(ax =ax)

plt.plot(pos[:,0], pos[:,1])
print pos.shape
plt.show()
