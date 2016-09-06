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
    def planet_plotter(self,positions, years, fig=None, ax = None):
        if not ax:
            fig, ax = plt.subplots()
        n = len(positions)
        print positions.shape
        hei = raw_input(" ")
        ax.plot(positions[:,0,:], positions[:,1,:])
        # (N, 7, 2)
        #colors = np.linspace(0,1,years+1)
        #ax.scatter(positions.T[0].T[::n/years], positions.T[1][::n/years])
        return fig, ax


    def simulate(self, years = 8, planet_index = "", res = 10000):
        G = 4*np.pi**2
        n = res*years
        stop = years
        dt = stop/float(n+1)

        starMass = self.starMass
        mass = self.mass

        times = np.linspace(0,stop,n)
        pos = np.zeros((n+1,2,self.numberOfPlanets))
        vel = np.zeros_like(pos)
        accel = np.zeros_like(pos)

        pos[0] = np.array((self.x0, self.y0))
        vel[0] = np.array((self.vx0, self.vy0))

        r = np.sqrt(pos[0,0]**2 +pos[0,1]**2)

        accel[0] = -G*starMass*pos[0]/r**3
        print "Simulating %d years in %d timesteps"%(years, n)

        for i in range(n): 
            if i%(n/100) == 0:
                print 100*i/float(n)
            r = np.sqrt(pos[i][0]**2 + pos[i][1]**2)
            accel[i] = -G*starMass*pos[i]/r**3
            vel[i+1] = vel[i] + accel[i] * dt
            pos[i+1] = pos[i] + vel[i+1] * dt
        return pos, accel

    def polar_orbits(self, planets= '', ax = None, plot = False, color=''):
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
        if color:
            ax.plot(x.T,y.T, c=color)
        else:
            ax.plot(x.T,y.T)
        #ax.plot(theta, orbits)

        if plot:
            plt.axis('equal')
            plt.show()
            return ax
        else:
            return ax



        
if __name__ == "__main__":
    ax = plt.subplot(111)
    seed =  87464 #adam:20776# fredrik:81995
    system = MySolarSystem(seed)
    ax = system.polar_orbits(ax=ax, color = '')
    plt.axis('equal')
    plt.show()

