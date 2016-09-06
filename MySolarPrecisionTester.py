import matplotlib.pyplot as plt
import numpy as np
from MySolarSystem import MySolarSystem

"""
File: MySolarPrecisionTester.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Supposed to calculate radius as a function of theta
numerically and analytically and show the difference. Not working
"""

def analytical(self, res = 250):
    indexes = range(self.numberOfPlanets)
        

    orbits = np.zeros((res, len(indexes)))
    theta = np.linspace(0,2*np.pi,res)
    for i in indexes:
        a = self.a[i]
        e = self.e[i]
        f = theta - self.psi[i]
        for j in range(res):
            orbits[j][i] = a*(1-e**2)/(1-e*np.cos(f[j]))

    return orbits, f

seed = 87464
system = MySolarSystem(seed)
years = 30; res = 250
pos, accel = system.simulate(years = years, res = res)


r = np.sqrt(pos[:,0]**2 +pos[:,1]**2)
theta = np.linspace(0,2*np.pi,res*years+1)
r2, theta2 =  analytical(system)

plt.plot(theta, r)
plt.plot(theta2, r2)

#ax.scatter(0,0, c='y')
#plt.axis('equal')
#plot_max = max(system.a +1)
#plt.axis([-plot_max,plot_max,-plot_max,plot_max])
#ax = system.polar_orbits(ax =ax)
plt.show()
