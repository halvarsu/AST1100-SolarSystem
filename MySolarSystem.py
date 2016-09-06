# coding: utf-8
from AST1100SolarSystem import AST1100SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class MySolarSystem(AST1100SolarSystem):

    """This is a doc string"""

    def __init__(self, seed):
        """TODO: to be defined1.

        :seed: TODO

        """
        AST1100SolarSystem.__init__(self, seed)

        self._seed = seed

    def polar_orbits(self, planets= '', ax = None, plot = False):
        '''Plots the orbits of the planets in AST1100SolarSystem given by
        the analytical solutions to the two body problem.

        Variable 'planets' should contain indexes of desirable planets
        seperated by commas. Empty string = all planets'''
        res = 250
        if not planets:
            indexes = range(self.numberOfPlanets)
        else:
            indexes = planets.split(',')
            

        orbits = np.zeros((res, len(indexes)))
        theta = np.linspace(0,2*np.pi,res)
        for i in indexes:
            a = self.a[i]
            e = self.e[i]
            f = theta - self.psi[i]
            for j in range(res):
                orbits[j][i] = a*(1-e**2)/(1-e*np.cos(f[j]))

        if not ax:
            ax = plt.subplot(111)#, projection = 'polar')

        x = np.array([r*np.cos(theta) for r in orbits.T])
        y = np.array([r*np.sin(theta) for r in orbits.T])

        ax.scatter(0,0, c='y')
        ax.plot(x.T,y.T)
        #ax.plot(theta, orbits)

        if plot:
            plt.axis('equal')
            plt.show()
            return ax
        else:
            return ax



        
if __name__ == "__main__":
    seed = 87464
    system = MySolarSystem(seed)
    system.polar_orbits(plot=True)

