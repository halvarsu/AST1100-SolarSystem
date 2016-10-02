from MySolarSystem import MySolarSystem 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import PyKEP as pkp
from PyKEP import epoch, AU, MU_SUN, DEG2RAD, lambert_problem, DAY2SEC, DAY2YEAR
from PyKEP.orbit_plots import plot_planet, plot_lambert
from MySateliteSim import MySateliteSim

class LambertSolarSystem(MySolarSystem):

    """A subclass of MySolarSystem where the module PyKEP is used. One of
    the uses is optimization of the MGA problem.
    """

    def __init__(self, seed):
        """TODO:

        :seed: The seed to be initialized in MySolarSystem

        """
        MySolarSystem.__init__(self, seed)

        self._seed = seed
        
    def init_planet(self, planet = 0, name = 'home', time = 0):
        '''Initializes an AST1100SolarSystem planet with the PyKEP module'''

        a,e,i,W,w,M = self.getKeplerianElements(new_values = True).T[planet]
        if time:
            M = self.getM(time, planet)
        mu_star = self.starMass * MU_SUN
        mu_planet = self.mass[planet] * MU_SUN 
        radius = self.radius[planet] * 1000
        safe_radius = radius + 100000
        return pkp.planet.keplerian(epoch(0), np.array((a*AU, e, i,
                W, w, M)), mu_star, mu_planet, radius, safe_radius,name)

    def lambert_solve(self, t1,t2,A,B, nameA='herbin', nameB='hars'):
        '''Solves the lambert problem in the AST1100SolarSystem system'''
        t1 = epoch(t1/DAY2YEAR)
        t2 = epoch(t2/DAY2YEAR)
        dt = (t2.mjd2000 - t1.mjd2000)*DAY2SEC

        planet1 = self.init_planet(A, name = nameA)
        planet2 = self.init_planet(B, name = nameB)

        r1, v1 = planet1.eph(t1)
        r2, v2 = planet2.eph(t2)
        #print np.array(r1) / AU

        mu_star = self.starMass * MU_SUN
        l = lambert_problem(r1, r2, dt, mu_star)
        return l

    def lambert_plot(self,a,A,B):
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        axis = fig.gca(projection = "3d")

        t1 = epoch(t1/DAY2YEAR)
        t2 = epoch(t2/DAY2YEAR)
        dt = (t2.mjd2000 - t1.mjd2000)*DAY2SEC

        planet1 = self.init_planet(self, A)
        plot_planet(planet1, t0=t1, legend=True, units = AU, ax = axis)
        planet4 = self.init_planet(self, B)
        plot_planet(planet4, t0=t2, legend=True, units = AU, ax = axis)

        rH, vH = planet1.eph(t1)
        rT, vT = planet4.eph(t2)

        mu_star = self.starMass * MU_SUN
        l = lambert_problem(rH, rT, dt, mu_star)
        plot_lambert(l, legend=True, units = AU, ax = axis)
        #print j
        #plot_lambert(l, sol=1, legend=True, units = AU, ax = axis)
        axis.set_xlabel('x')
        plt.show()
        return l
