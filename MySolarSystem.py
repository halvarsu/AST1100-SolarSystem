# coding: utf-8
from AST1100SolarSystem import AST1100SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

ang2pix = AST1100SolarSystem.ang2pix

class MySolarSystem(AST1100SolarSystem):
    """A subclass of AST1100SolarSystem made for easier implementation of
    own methods using data from the previous mentioned module. """

    def __init__(self, seed):
        """TODO: 

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


    def saveData(self, times, pos, vel=None, accel=None):
        np.save('data/timesData', times)
        np.save('data/positionsData', pos)
        np.save('data/velocityData', vel)
        np.save('data/accelerationData', accel)

    def analytical_pos(self, planets= '', res = 1000, xy = False):
        # If no index supplied, uses all planets
        if not planets:
            indexes = range(self.numberOfPlanets)
        # If index is int, check if valid int
        elif type(planets) == int and planets <= self.numberOfPlanets:
            indexes = str(planets)
        else:
            indexes = planets.split(',')
            indexes = np.array([int(i) for i in indexes])

        orbits = np.zeros((res, len(indexes)))
        theta = np.linspace(0,2*np.pi,res)
        a_arr = self.a[indexes]
        e_arr = self.e[indexes]
        for i in range(len(indexes)):
            a = a_arr[i]
            e = e_arr[i]
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

    def getExactPos(self):
        if not hasattr(self, 'exactPos'):
            exactPos = np.load('positionsHomePlanet.npy')
            self.exactPos = exactPos
            return self.exactPos
        return self.exactPos

    def getTimes(self):
        if not hasattr(self,'times'):
            times = np.load('times.npy')
            self.times = times
            return self.times
        return self.times

    def getTheta(self):
        if not hasattr(self, 'theta') or not hasattr(self,'orbits_r'):
            self.theta, self.orbits_r = self.analytical_pos()
        return self.theta, self.orbits_r

    def getPositionsFunction(self, planets = '', xy = False):
        '''  posFunction returns positions as function of time.'''
        if not hasattr(self, 'posFunction'):
            from scipy.interpolate import interp1d
            #from scipy.interpolate import UnivariateSpline
            planetPos   = self.getExactPos()
            times    = self.getTimes()
            self.posFunction = interp1d(times[:-1], planetPos)
            #self.posFunction = UnivariateSpline(times[:-1], planetPos)
        return self.posFunction 

    def getAngleFunction(self):
        '''angleFunction returns positions as function of theta'''
        if not hasattr(self, 'angleFunction'):
            from scipy.interpolate import interp1d
            theta, orbits_r = self.getTheta()
            angleFunction = interp1d(theta, orbits_r.T)
        return angleFunction

    def velFunction(self,t,dt = 0.000001):
        return (self.posFunction(t+dt) - self.posFunction(t))/dt

    def getVelocityFunction(self):
        'velfunction needs posfunction to function'
        if not hasattr(self, 'posFunction'):
            self.getPositionsFunction()
        if not hasattr(self, 'velFunction'):
            print "Something is wrong, velFunction not implemented!"
        return self.velFunction

    def getKeplerianElements(self, new_values = False):
        '''Returns the six orbital elements needed by PyKEP for
        instantiation of planets. The elements are
        a : Semi major axis
        e : Eccentricity
        i : Inclination  ->  0  (2-dim)
        W : Longitude of ascending node -> 0 (2-dim)
        w : Argument of periapsis
        M : Mean anomaly 
        Returned array has shape (6, numberOfPlanets)'''

        a = self.a  
        e = self.e
        inc = np.zeros(self.numberOfPlanets)
        W = np.zeros(self.numberOfPlanets)
        try:
            if hasattr(self, 'w'):
                w = self.w
            elif new_values:
                raise IOError
            else:
                w = np.load('data/orbitalData/argumentOfPeriapsis.npy')
        except IOError:
            print "Creating new values"
            pos = self.getExactPos()
            r = np.linalg.norm(pos, axis = 0)
            index_periapsis = np.argmin(r, axis = 1)
            pos_periapsis   = pos[:,:,index_periapsis]
            b = a*np.sqrt(1-e**2)
            w = np.zeros(self.numberOfPlanets)
            for i in range(len(w)):
                print a[i], b[i], pos_periapsis[0,i,i], pos_periapsis[1,i,i]
                w[i] =\
                np.arctan((a[i]*pos_periapsis[1,i,i])/(b[i]*pos_periapsis[0,i,i]))
                # actan only gives values for [-pi/2 to pi/2]
                w[i] += np.pi if pos_periapsis[0,i,i] < 0 else 0
            np.save('data/orbitalData/argumentOfPeriapsis.npy', w)
        try:
            if hasattr(self, 'w'):
                M = self.M
            elif new_values:
                raise IOError
            else:
                M = np.load('data/orbitalData/meanAnomaly.npy')
        except IOError:
            b = a*np.sqrt(1-e**2)
            x0 = self.x0; y0 = self.y0
            M = np.arctan((y0*a)/(x0*b)) - e*y0/b 
            M -= w
            M[x0<0] += np.pi
            np.save('data/orbitalData/meanAnomaly.npy', M)

        self.M = M
        self.w = w
        data = np.array((a,e,inc,W,w,M))
        return data

    def getM(self, t):
        posFunc = self.getPositionsFunction()
        a = self.a; e = self.e
        b = a*np.sqrt(1-e**2)
        x, y = posFunc(t)
        M = np.arctan((y*a)/(x*b)) - e*y/b 
        M -= self.w
        M[x<0] += np.pi
        return M





    def find_hohmann_params(self, A=0, B=1):
        '''returns parameters of the hohmann transfer from planet A to planet B
        for an approximation of planet orbits as circles
        
        : A, B : : Expects B > A
        '''
        pi = np.pi
        posFunction, angleFunction = self.find_functions("%d, %d" %(A,B))
        n = 10
        r = np.zeros((n,2))
        for i,thet in enumerate(np.linspace(0,2*pi, n)):
            r[i] = angleFunction(thet)

        rA = np.mean(r, axis = 0)[0]
        rB = np.mean(r, axis = 0)[1]

        aTransfer = (rA+rB)/2.
        G = 4*pi**2
        M = self.starMass
        v_init_A        = np.sqrt(G*M/rA)
        v_final_B       = np.sqrt(G*M/rB)
        v_transfer_A    = np.sqrt(G*M*(2./rA - 1./aTransfer))
        v_transfer_B    = np.sqrt(G*M*(2./rB - 1./aTransfer))
        delta_VA        = v_transfer_A - v_init_A
        delta_VB        = v_final_B - v_transfer_B
        total_delta_V   = delta_VA + delta_VB
        e               = 1 -     rA/aTransfer #eccentrisity of transfer ellipse
        transferTime    = pi*np.sqrt((rA+rB)**3/(8*G*M))

        return {"v_init_A": v_init_A, 
                "v_final_B": v_final_B, 
                "total_delta_V":total_delta_V, 
                "e":e, 
                "transferTime":transferTime,
                "v_transfer_A":v_transfer_A,
                "v_transfer_B":v_transfer_B}
        


        
if __name__ == "__main__":
    ax = plt.subplot(111)
    seed =  87464 #adam:20776# fredrik:81995
    system = MySolarSystem(seed)
    ax = system.polar_orbits(ax=ax, color = '')
    plt.axis('equal')
    plt.show()

