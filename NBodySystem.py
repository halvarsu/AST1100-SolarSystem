from MySolarSystem import MySolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class NBodySystem(MySolarSystem):

    """Docstring for NBodySystem. """

    def __init__(self, seed):
        """TODO: to be defined1.

        :seed: TODO

        """
        MySolarSystem.__init__(self, seed)
        self.data_filenames = 'data/{}NBody200yr{}'

        self._seed = seed
            
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
        return pos, vel, accel, times

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

    def loadData(self):
        self.pos = np.load(self.data_filenames.format('pos','.npy'))
        self.vel = np.load(self.data_filenames.format('vel','.npy'))
        self.acc = np.load(self.data_filenames.format('acc','.npy'))
        self.times = np.load(self.data_filenames.format('times','.npy'))
        return self.pos, self.vel, self.acc, self.times


    def saveData(self, times, pos, vel=None, acc=None):
        np.save(self.data_filenames.format('pos',''), pos)
        np.save(self.data_filenames.format('acc',''), acc)
        np.save(self.data_filenames.format('vel',''), vel)
        np.save(self.data_filenames.format('times',''), times)


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
