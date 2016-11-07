# coding: utf-8
from AST1100SolarSystem import AST1100SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_engine_data(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    data = [int(lines[0])] # first line, outside_count
    for line in lines[1:]:
        data.append(float(line))
    return data

ang2pix = AST1100SolarSystem.ang2pix

class MySolarSystem(AST1100SolarSystem):
    """A subclass of AST1100SolarSystem made for easier implementation of
    own methods using data from the previous mentioned module. """

    au = 149597870700 # m
    year = int(365.25*24*3600) #s
    day2sec = float(24*3600) #s

    def __init__(self, seed):
        """TODO: Rework documentation 
        :seed: The seed to be initialized in AST1100SolarSystem
        """
        AST1100SolarSystem.__init__(self, seed)

        self._seed = seed
    def simulate(self, years = 8, planet_index = "", res = 10000,
            integrator = 'leapfrog'):
        '''
        Simulates the orbits of the planets given by python module
        AST1100SolarSystem. 

        Parameters
        ----------
        years : float
            Number of years to run simulation

        planet_index : str with comma seperated ints
            String containing indexes of planets uses.

        res : int
            Amount of timesteps per year

        integrator : str
            Integrator to use. Either 'leapfrog' or 'euler cromer'
        '''
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
            print "JEA"
            for i in xrange(n):
                if i%(n/100) == 0:
                    print 100*i/float(n)
                r = np.linalg.norm(pos[i],axis=0)
                pos[i+1] = pos[i] + vel[i]*dt + 0.5*accel[i]*dt**2
                accel[i+1] = -G*starMass*pos[i+1]/r**3
                vel[i+1] = vel[i] + 0.5*(accel[i]+accel[i+1])*dt
        elif integrator == 'euler_cromer':
            for i in range(n): 
                if i%(n/100) == 0:
                    print 100*i/float(n)
                r = np.linalg.norm(pos[i],axis=0)
                accel[i] = -G*starMass*pos[i]/r**3
                vel[i+1] = vel[i] + accel[i] * dt
                pos[i+1] = pos[i] + vel[i+1] * dt
        else:
            print "integrator %s not implemented" %integrator
            assert NotImplementedError, "not implemented"
            import sys
            sys.exit(1)
        self.saveData(times, pos, vel, accel)
        return pos, vel, accel, times


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

        Parameters
        ----------

        planets : str with comma seperated values
            planets to find orbits of. empty string = all planets

        ax : matplotlib axes
            Axes to plot the orbits onto.

        plot : bool
            Whether to plot or not. Equalizes axes.

        color : str
            color of the orbits
        '''

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
        """returns self.exactPos, loaded from positionsHomePlanet.npy"""
        if not hasattr(self, 'exactPos'):
            exactPos = np.load('positionsHomePlanet.npy')
            self.exactPos = exactPos
            return self.exactPos
        return self.exactPos

    def getTimes(self):
        """returns self.times, loaded from times.npy"""
        if not hasattr(self,'times'):
            times = np.load('times.npy')
            self.times = times
            return self.times
        return self.times
        folder = 'data/atmosphericData/density'
        folder = 'data/atmosphericData/density'

    def getTheta(self):
        if not hasattr(self, 'theta') or not hasattr(self,'orbits_r'):
            self.theta, self.orbits_r = self.analytical_pos()
        return self.theta, self.orbits_r

    def getDensityFunction(self):
        if not hasattr(self, 'angleFunction'):
            from scipy.interpolate import interp1d
            folder = 'data/atmosphericData/density'
            h = np.load(folder+'height.npy')
            rho = np.load(folder+'density.npy')
            self.densityFunc = interp1d(h, rho, bounds_error = False, 
                                        fill_value = 0)
        return self.densityFunc

    def getPressureFunction(self):
        if not hasattr(self, 'angleFunction'):
            from scipy.interpolate import interp1d
            folder = 'data/atmosphericData/density'
            h = np.load(folder+'height.npy')
            P = np.load(folder+'pressure.npy')
            self.pressureFunc = interp1d(h, P, bounds_error = False, 
                                         fill_value = 0)
        return self.pressureFunc

    def getTemperatureFunction(self):
        if not hasattr(self, 'angleFunction'):
            from scipy.interpolate import interp1d
            folder = 'data/atmosphericData/density'
            h = np.load(folder+'height.npy')
            T = np.load(folder+'temperature.npy')
            self.tempFunc = interp1d(h, T, bounds_error = False, 
                                     fill_value = 130)
        return self.tempFunc


    def getPositionsFunction(self, planets = '', xyz = False, 
                             force_new=False):
        '''Loads and returns posFunction, which gives position of the
        planets as function of time.'''
        if not hasattr(self, 'posFunction') or force_new:
            from scipy.interpolate import interp1d
            #from scipy.interpolate import UnivariateSpline
            if xyz:
                planetPos2d = self.getExactPos()
                s = planetPos2d.shape
                planetPos = np.zeros((3,s[1],s[2]))
                planetPos[:-1,:,:] = planetPos2d
            else:
                planetPos = self.getExactPos()
            times = self.getTimes()
            self.posFunction = interp1d(times[:-1], planetPos)
            #self.posFunction = UnivariateSpline(times[:-1], planetPos)
        return self.posFunction 

    def getAngleFunction(self):
        '''angleFunction returns positions of the planets as function 
        of theta'''
        if not hasattr(self, 'angleFunction'):
            from scipy.interpolate import interp1d
            theta, orbits_r = self.getTheta()
            angleFunction = interp1d(theta, orbits_r.T)
        return angleFunction

    def velFunction(self,t,dt = 0.000001):
        return (self.posFunction(t+dt) - self.posFunction(t))/dt

    def getVelocityFunction(self):
        '''velfunction needs posfunction to function'''
        if not hasattr(self, 'posFunction'):
            self.getPositionsFunction()
        if not hasattr(self, 'velFunction'):
            print "Something is wrong, velFunction not implemented!"
        return self.velFunction

    def testVelFunctionAccuracy(self):
        velFunc= self.getVelocityFunction()
        exact = np.linalg.norm((self.vx0[0], self.vy0[0]))
        for i in range(10):
            test = np.linalg.norm(velFunc(0, 10**-i), axis = 0)[0]
            print test, abs((test-exact)/exact), "   ", 10**-i
        print exact

    def getKeplerianElements(self, new_values = False):
        '''Returns the six orbital elements needed by PyKEP for
        instantiation of planets. The elements are
        a : Semi major axis
        e : Eccentricity
        i : Inclination  ->  0  (2-dim)
        W : Longitude of ascending node -> 0 (2-dim)
        w : Argument of periapsis
        M : Mean anomaly 
        Returned array has shape (6, numberOfPlanets).
        DISCLAIMER: PyKEP has not produces the right
        starting points for with data from this method as of today. 
        Reason unknown. '''

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
                np.arctan( (a[i]*pos_periapsis[1,i,i])/(
                    b[i]*pos_periapsis[0,i,i]))
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

    def initConstants(self):
        """
        Initializes 
        -----------
        self.c : float
            the speed of light 

        self.h : float
            Plancks constant 

        self.k_B : float
            Boltzmans constant 
        
        Returns
        -------
        c, h, k_B
        """
        if not hasattr(self, 'k_B'):
            import astropy.constants as const
            self.c = const.c.value      # Speed of light 
            self.h = const.h.value      # Plancks constant?
            self.k_B = const.k_B.value  # Boltzmans constant?
        return self.c, self.h, self.k_B

    def wavelengthIntensity(self, wl, T = None):
        """Returns intensity as function of wavelength and temp for a 
        black body
        Parameters
        ----------
        wl : float or sequence of floats
            Wavelength
        T : float
            Temperature"""
        if not T:
            T = self.temperature
        h, c, k_B = self.initConstants()
        return 2*h*c**2/wl**5 * 1/(np.exp(h*c/(k_B*T*wl))-1)

    def plotBlackbody(self, T = None):
        """Assumes a black body star, plots intensity as function of
        wavelengths and the given temperature.
        Parameters
        ----------
        T : float or sequence of floats
            Temperature of black body

        Returns
        -------
        out : flux. 
            Flux of black body. Not tested"""
        if not T:
            T = self.temperature
        self.initConstants()
        wavelengths = np.linspace(1,3000,200)
        intensity = self.wavelengthIntensity(wavelengths*1e-9, T)

        plt.plot(wavelengths, intensity)
        plt.title('T = %d K' %int(T))
        plt.xlabel('Wavelength $\lambda$ in nm')
        plt.ylabel('Intensity B($\lambda$)')
        plt.savefig('figure/wavelength_plot.png')
        plt.show()

        self.flux = 2*np.pi**5*k**4/(15*h**3*c**2)*T**4
        print " Flux = %g"% self.flux # not tested properly
        return self.flux

    def find_hohmann_params(self, A=0, B=1):
        ''' WORK IN PROGRESS
        returns parameters of the hohmann transfer from planet A 
        to planet B for an approximation of planet orbits as circles.

        Parameters
        ----------
        A, B : int
            Index of planets A and B. Expects B > A
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

    def get360Projections(self):
        self.proj = np.load('all_the_projections.npy')
        return self.proj

        
    def projection(self, phi_0, theta_0=np.pi/2):
        from PIL import Image
        from numpy import cos, sin, pi, arcsin, arctan, sqrt
        pi = np.pi; cos = np.cos; sin = np.sin
        arcsin = np.arcsin; arctan = np.arctan; sqrt = np.sqrt 
        ang2pix = AST1100SolarSystem.ang2pix

        folder = 'data/projections/'
        inFile = open(folder+'himmelkule.npy', 'rb')
        celestial_sky = np.load(inFile)
        inFile.close()

        height = 480
        width  = 640 
        fov_p = 2*pi*70/360.    #radians
        fov_t = 2*pi*70/360.    #radians
        xlim = 2*sin(fov_p/2.)/(1+cos(fov_p/2.)) 
        ylim = 2*sin(fov_t/2.)/(1+cos(fov_t/2.)) 
        x_pic = np.linspace(-xlim, xlim, width)
        y_pic = np.linspace(ylim, -ylim, height)
        proj_rgb = np.zeros((height,width,3), dtype='uint8')

        for i,x in enumerate(x_pic):
            for j,y in enumerate(y_pic):
                rho = sqrt(x**2 + y**2)
                c = 2 * arctan(rho/2)
                if rho == 0:
                    phi = phi_0
                    theta = theta_0
                else:
                    theta = pi/2-arcsin(
                                    cos(c)*cos(theta_0) + \
                                    y*sin(c)*sin(theta_0)/rho)
                    phi = phi_0 + arctan(
                                    x*sin(c)/(rho*sin(theta_0)*cos(c) -\
                                    y*cos(theta_0)*sin(c)))
                pix = ang2pix(theta, phi)
                proj_rgb[j,i,:] = celestial_sky[pix][2:]
                    
        return proj_rgb

    def getPlanetTemp(self):
        T_star = self.temperature
        r = self.a * self.au
        R_star = self.starRadius * 1000
        self.planetTemperature = T_star*np.sqrt(R_star/(np.sqrt(2)*r))
        return self.planetTemperature


    def pressure_solver(self):
        from scipy.constants import m_p, m_e, k, G
        import matplotlib.pyplot as plt

        sunMass = 1.989e30
        mu = 38
        m_H = m_p + m_e
        r_p = self.radius[5] * 1000
        M = self.mass[5] * sunMass

        rho0 = self.rho0[5] 
        T0 = self.getPlanetTemp()[5]
        print T0
        P0 = rho0 * k * T0 / (mu * m_H)
        gamma = 1.4
        beta = P0**(1-gamma)*T0**gamma
        n = 300000

        r = np.linspace(r_p,r_p+250000,n)
        P = np.zeros(n)
        rho = np.zeros_like(P)
        T = np.zeros_like(P)
        P[0] = P0
        rho[0] = rho0
        T[0] = T0

        dr = r[1] - r[0]
        i= 0
        tol = 1e-7
        while P[i] > P0*tol:
            if i >= (n-1):
                print "Breaking early"
                break
            P[i+1] = P[i] - rho[i]*G*M/r[i]**2*dr
            if T[i] > T0/2:
                T[i+1] = beta**(1./gamma)/P[i+1]**((1-gamma)/gamma)
            else:
                T[i+1] = T[i]
            rho[i+1] = P[i+1]*mu*m_H/(k*T[i+1])
            i += 1

        plt.plot(r[:i]-r_p, rho[:i])
        plt.title('Density')
        plt.xlabel('Height [m]')
        plt.ylabel('density $\\rho$ [kg/m^3]')

        plt.show()
        folder = 'data/atmosphericData/density'
        np.save(folder+'height', r[:i] - r_p)
        np.save(folder+'density', rho[:i])
        np.save(folder+'pressure', P[:i])
        np.save(folder+'temperature', T[:i])




if __name__ == "__main__":
    ax = plt.subplot(111)
    seed =  87464 #adam:20776# fredrik:81995
    system = MySolarSystem(seed)
    ax = system.polar_orbits(ax=ax, color = '')
    plt.axis('equal')
    plt.show()

