# coding: utf-8
from AST1100SolarSystem import AST1100SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class MySolarSystem(AST1100SolarSystem):

    """A subclass of AST1100SolarSystem made for easier implementation of
    own methods using data from the previous mentioned module. """

    def __init__(self, seed):
        """TODO: to be defined1.

        :seed: The seed to be initialized in AST1100SolarSystem

        """
        AST1100SolarSystem.__init__(self, seed)

        self._seed = seed
    def planet_plotter(self,positions, years, fig=None, ax = None):
        if not ax:
            fig, ax = plt.subplots()
        n = len(positions)
        print positions.shape
        ax.plot(positions[:,0,:], positions[:,1,:])
        # (N, 7, 2)
        colors = np.linspace(0,1,years+1)
        ax.scatter(positions[:,0,:][::n/years],
                positions[:,1,:][::n/years])#, c = colors)
        return fig, ax

    def simulate(self, years = 8, planet_index = "", res = 10000,
            integrator = 'leapfrog'):
        G = 4*np.pi**2
        n = res*years
        stop = years
        dt = stop/float(n+1)

        starMass = self.starMass
        mass = self.mass

        times = np.linspace(0,stop,n+1)
        pos = np.zeros((n+1,2,self.numberOfPlanets))
        vel = np.zeros_like(pos)
        accel = np.zeros_like(pos)

        pos[0] = np.array((self.x0, self.y0))
        vel[0] = np.array((self.vx0, self.vy0))

        r = np.sqrt(pos[0,0]**2 +pos[0,1]**2)

        accel[0] = -G*starMass*pos[0]/r**3
        print "Simulating %d years in %d timesteps"%(years, n)

        if integrator == 'leapfrog':
            for i in xrange(n):
                if i%(n/100) == 0:
                    print 100*i/float(n)
                r = np.sqrt(pos[i][0]**2 + pos[i][1]**2)
                pos[i+1] = pos[i] + vel[i]*dt + accel[i]*dt**2
                accel[i+1] = -G*starMass*pos[i]/r**3
                vel[i+1] = vel[i] + 0.5*(accel[i]+accel[i+1])*dt

        elif integrator == 'euler_cromer':
            for i in range(n): 
                if i%(n/100) == 0:
                    print 100*i/float(n)
                r = np.sqrt(pos[i][0]**2 + pos[i][1]**2)
                accel[i] = -G*starMass*pos[i]/r**3
                vel[i+1] = vel[i] + accel[i] * dt
                pos[i+1] = pos[i] + vel[i+1] * dt
        else:
            print "integrator %s not implemented" %integrator
            assert NotImplementedError, "not implemented"
            import sys
            sys.exit(1)
        self.saveData(times, pos, vel, accel)
        return pos, accel, times

    def simulate_nbody(self, years = 8, planet_index = "", res = 5000):
        self.G = 4*np.pi**2
        n = res*years
        numOfBodies = self.numberOfPlanets + 1
        stop = years
        dt = stop/float(n+1)


        times = np.linspace(0,stop,n+1)
        pos = np.zeros((n+1,2,numOfBodies))
        vel = np.zeros_like(pos)
        accel = np.zeros_like(pos)

        self.star_velocity_init()
        mass_array = np.append(self.mass,self.starMass)
        x0 = np.append(self.x0, self.star_x0)
        y0 = np.append(self.y0, self.star_y0)
        vx0 = np.append(self.vx0, self.star_vx0)
        vy0 = np.append(self.vy0, self.star_vy0)
        
        pos[0] = np.array((x0, y0))
        vel[0] = np.array((vx0, vy0))

        accel[0] = self.accel_nbody(pos[0], mass_array)
        print "Simulating %d years in %d timesteps"%(years, n)

        for i in xrange(n):
            if i%(n/100) == 0:
                print 100*i/float(n)
            pos[i+1] = pos[i] + vel[i]*dt + accel[i]*dt**2
            accel[i+1] = self.accel_nbody(pos[i+1], mass_array)
            vel[i+1] = vel[i] + 0.5*(accel[i]+accel[i+1])*dt
        self.saveData(times, pos, vel, accel)
        return pos, accel, times

    def star_velocity_init(self):
        """ 
        Sets the velocity of the star so it matches the movement of the
        planets (ensures a stationary center of mass for nbody-sim).
        """
        starMass = self.starMass
        masses = self.mass
        self.star_vx0 = 0
        self.star_vy0 = 0
        self.star_x0 = 0
        self.star_y0 = 0

        for i in range(self.numberOfPlanets):
            self.star_vx0 -= self.vx0[i]*masses[i]/starMass
            self.star_vy0 -= self.vy0[i]*masses[i]/starMass
            self.star_x0 -=  self.x0[i]*masses[i]/starMass
            self.star_y0 -=  self.y0[i]*masses[i]/starMass


    def accel_nbody(self,pos,mass):
        """Used in simulate_nbody to find the acceleration of all the
        bodies on eachother"""
        N = len(mass)
        f_array = np.zeros((N,2))
        for i in range(N-1):
            for j in range(i+1,N):
                r_vec = pos[:,i] - pos[:,j]
                r = np.linalg.norm(r_vec)
                f = - self.G * mass[i]*mass[j] *r_vec/r**3
                f_array[i] += f
                f_array[j] -= f
        return (f_array.T/mass)
     

    def saveData(self, times, pos, vel, accel):
        np.save('data/timesData', times)
        np.save('data/positionsData', pos)
        np.save('data/velocityData', vel)
        np.save('data/accelerationData', accel)

    def analytical_pos(self, planets= '', res = 1000, xy = False):
        if not planets:
            indexes = range(self.numberOfPlanets)
        else:
            indexes = planets.split(',')
            indexes = [int(i) for i in indexes]
            
        print 'Finding analytical solution for these planets:',indexes

        orbits = np.zeros((res, len(indexes)))
        theta = np.linspace(0,2*np.pi,res)
        for i in indexes:
            a = self.a[i]
            e = self.e[i]
            f = theta - self.psi[i]
            for j in range(res):
                orbits[j][i] = a*(1-e**2)/(1-e*np.cos(f[j]))
        if xy:
            x = np.array([r*np.cos(theta) for r in orbits.T])
            y = np.array([r*np.sin(theta) for r in orbits.T])
            return x, y
        else: 
            return theta, orbits


    def polar_orbits(self, planets= '', ax = None, plot = False, color=''):
        '''Plots the orbits of the planets in AST1100SolarSystem given by
        the analytical solutions to the two body problem.

        Variable 'planets' should contain indexes of desirable planets
        seperated by commas. Empty string = all planets'''
        res = 250
        x,y= self.analytical_pos(planets, res = res, xy = True)
        if not ax:
            ax = plt.subplot(111)#, projection = 'polar')

        ax.scatter(0,0, c='y')
        if color:
            ax.plot(x.T,y.T, c=color)
        else:
            #ax.scatter(x.T[0], y.T[0])
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

