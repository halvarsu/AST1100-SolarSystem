# coding: utf-8
from MySolarSystem import MySolarSystem
import numpy as np
import matplotlib.pyplot as plt


ang2pix = MySolarSystem.ang2pix


def rot_matrix(angle):
    """Provides a two dimensional rotation matrix for a given angle"""
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array(((c,-s),(s,c)))

def find_distance(distances, time = 0):
    '''
    Finds satelite position from the distance to the planets at a 
    given time.
    Parameters
    ----------
    distances : list of floats, len x+1
        Distance to the x first planets and the star. Star should be 
        last element
    '''

    r_s = distances[-1]
    rp_list = distances[:-1]
    n = len(rp_list) 
    system = m.MySateliteSim(87464)

    posFunc = system.getPositionsFunction()
    planetPos = posFunc(time).T
    print planetPos.shape

    theta = np.linspace(0,2*np.pi, 10000)
    circle = r_s* np.array((np.sin(theta), np.cos(theta)))
    plt.plot(circle[0], circle[1],color = 'k')
    for i in range(n):
        x0 = planetPos[i,0]
        y0 = planetPos[i,1]
        rp = rp_list[i]
        print '---------------------------'

        A = (x0**2 + y0**2 + r_s**2 - rp**2)/(2*x0)
        b = y0/float(x0)
        C = np.sqrt(r_s**2*(b**2+1) - A**2)
        y = [(A*b + C)/(b**2+1),(A*b - C)/(b**2+1)]
        x = [A- y_i*b for y_i in y]
        for i in range(2):
            plt.scatter(x[i], y[i])
            print "comp", x[i], y[i] 
        circle = np.array((rp*np.sin(theta)+x0, rp*np.cos(theta)+y0))
        plt.plot(circle[0], circle[1]) 
    plt.scatter(0,0,c='y', s=100)
    plt.axis('equal')
    plt.show()


def get_sat_velocity(dl_obs1, dl_obs2):
    '''Receives the doppler shift of the two reference stars and 
    returns the velocity of the spacecraft in the reference system 
    of the home star.
    Parameters
    ----------
    dl_obs1 : float
        Observed doppler shift of ref star 1 in m
    dl_obs2 : float
        Observed doppler shift of ref star 2 in m
    Returns
    -------
    v1, v2 : floats
        velocity of spacecraft, m/s '''

    c = const.c.value
    h_alpha = 656.3e-9
    dl_ref1 = 0.020225790030e-9
    phi1 = 190.439790
    dl_ref2 = -0.013132758138e-9
    phi2 =  118.096595
    v_refstar1 = c * dl_ref1 / h_alpha 
    v_refstar2 = c * dl_ref2 / h_alpha
    v1_obs = c * dl_obs1 / h_alpha 
    v2_obs = c * dl_obs2 / h_alpha 
    v_sat1 = v_refstar1 - v1_obs 
    v_sat2 = v_refstar2 - v2_obs

    phi1_rad = phi1 / 360. *2*np.pi
    phi2_rad = phi2 / 360. *2*np.pi
    s1 = np.sin(phi1_rad);  c1 = np.cos(phi1_rad)
    s2 = np.sin(phi2_rad);  c2 = np.cos(phi2_rad)
    inv_rel = 1/(np.sin(phi2_rad-phi1_rad))*np.array(((s2, -s1),(-c2,c1)))
    v1, v2 = np.dot(inv_rel, np.array((v_sat1, v_sat2)))
    return v1, v2



class MySateliteSim(MySolarSystem):
    """This is an implementation of My Solar System designed for sending
    satelites into space"""
    au = 149597870700 # m
    year = int(365.25*24*3600) #s

    def __init__(self, seed): 
        """
        On init, creates an instance of MySolarSystem with the seed
        provided

        :seed: The seed to be instantiated in MySolarSystem

        """
        MySolarSystem.__init__(self, seed)

        self._seed = seed
        self.initiated  =  False


    def initFunctions(self):
        """ loads arrays with planet positions and times.  These are then 
        used to create interpolation functions for position, velocity and
        acceleration of the planets.  """
        self.planetPos  =  np.load('positionsHomePlanet.npy')
        self.times      =  np.load('times.npy')
        self.posFunc    =  self.getPositionsFunction()
        self.velFunc    =  self.getVelocityFunction()
        self.angleFunc  =  self.getAngleFunction()
        self.accelFunc  =  self.getAccelFunction()


    def accelFunction(self, pos, time):
        G = 4*np.pi**2
        mass = self.mass
        r = np.sqrt(pos[0]**2 + pos[1]**2)
        acc = -G*self.starMass*pos/r**3

        planetPositions = self.posFunc(time).T
        for i,posPlanet in enumerate(planetPositions):
            rel_pos = pos - posPlanet
            r2 = rel_pos[0]**2 + rel_pos[1]**2
            acc += -G*mass[i]*rel_pos/(r2*np.sqrt(r2))
        return acc

    def getAccelFunction(self):
        #if not hasattr(self, 'accelFunction'):
            #self.accelFunction = np.vectorize(self.accelFunc_unvec)
        return self.accelFunction


    def orient(self, filename):
        from PIL import Image
        in_image = Image.open(filename)
        proj360 = self.get360Projections()
        N = proj360.shape[0]
        sums = np.zeros(N)
        for i,proj in enumerate(proj360):
            sums[i] = np.sum(np.sum(np.sum((proj - in_image)**2)))
        degrees = np.arange(N)

        plt.plot(degrees, sums)
        plt.show()
        best_fit = np.argmin(sums) 
        return best_fit

    def test_orient(self):
        from PIL import Image
        exact = 2*np.pi*np.random.random()
        exact_deg = exact*180./np.pi
        proj = self.projection(exact)
        img = Image.fromarray(proj)
        test_name = 'orient_test.png'
        img.save(test_name)
        computed = self.orient(test_name)
        tol = 0.5
        if np.abs(exact_deg-computed) < tol:
            success = True
        elif np.abs(exact_deg-computed-360) < tol:
            success = True
        elif np.abs(exact_deg-computed+360) < tol:
            success = True
        else:
            success = False
        print "Exact: ", exact_deg, "computed: ", computed
        assert success, "mismatch between calculated and exact degrees"


    def boost(self, vel, boost, boost_angle=0, boost_dir = None):
        """
        Applies an instantanious boost at an angle boost_angle at the 
        vel direction
        Parameters
        ----------
        vel : array of two floats           
            Two dimensional velocity before boost 
        boost : float or int        
            Magnitude of boost
        boost_dir : array of two floats
            boost dir, optional. If None, dir of vel is used
        boost_angle : angle
            angle of boost from boost_dir
        """

        direction = vel/np.linalg.norm(vel)
        if hasattr(boost_dir, '__len__'):
            boost_direction = np.dot( rot_matrix(boost_angle), 
                                      np.array(boost_dir) )
        else:
            boost_direction = np.dot( rot_matrix(boost_angle), 
                                      direction )
        new_vel = vel + boost_direction*boost/float(self.au)*self.year
        return new_vel


    def boost_init(
            self, pos0, vel0, t0, init_boost = 0, boost_angle = 0,
            boost_dir = None):
        """Sets up satelite position pos0 and velocity vel0 at time t0 
        with the given boost applied. 
        """
        self.initiated = True
        self.start = t0
        self.vel0 = self.boost(vel0, init_boost, boost_angle, boost_dir)
        self.pos0 = pos0

    def launch_init(
            self, t0 = 14.32, init_boost = 9594, init_angle = 1.6*np.pi, 
            start_dist=1000, start_planet = 0, 
            ground_launch = False):

        """Sets up satelite position relative to another planet"""
        self.initiated = True
        self.start = t0
        self.initFunctions()
        planetVel0 = self.velFunc(t0)
        planetPos0 = self.posFunc(t0)
        planet_angle = np.arctan(
                planetPos0[1,start_planet]/planetPos0[0,start_planet])
        if planetPos0[0,start_planet] < 0:
            planet_angle  += np.pi
        if ground_launch:
            sp_vel = planetVel0[:,start_planet]
            planet_direction = sp_vel/np.linalg.norm(sp_vel)
            # init_angle should be 0 for a straight launch trajectory:
            theta = init_angle + planet_angle + np.pi/2 
            init_pos_km =  planet_direction * self.radius[start_planet]
            init_vel = planet_direction *init_boost
        else:
            # init_angle important for direction of hyperbolic excess vel
            theta = init_angle + planet_angle 
            init_vel = np.array((-init_boost*np.sin(theta),
                                  init_boost*np.cos(theta)))
            init_dist = start_dist+self.radius[start_planet]
            init_pos_km = np.array((init_dist*np.cos(theta),
                                    init_dist*np.sin(theta)))
        au = float(self.au)
        pos0 = planetPos0.T[start_planet] + init_pos_km/au*1000
        vel0 = planetVel0.T[start_planet] + init_vel*self.year/au
        self.pos0 = pos0
        self.vel0 = vel0

    def write_launch_data(self):
        if not self.initiated:
            raise NotImplementedError(
                    'no init method called before load')

    def satelite_sim(
        self, tN = 4, dt_ip = 1/(50000.), dt_close = 1/(365.25*24*3600), 
        speed_factor = [1,1,1], target_planet = 4, break_close = True, 
        *args, **kwargs):

        ''' 
        Simulates satelite trajectory . There are two options, either a
        first_launch from planet, or a continuation from a point in space.
        Parameters
        ----------
        tN : float
            Number of years to run simulation
        dt_ip : float
            Size of interplanetary timesteps
        dt_close : float
            Size of near planet timesteps
        speed_factor : list of floats, len = 3
            Multiplication factor for the timestep of the  three 
            standard legs of a journey: close to departure planet, between
            planets, close to arrival planet.
        target_planet : int
            index of the planet which will trigger break at closest
            approach
        break_close : bool 
            whether to break at closest approach
        ''' 
        
        if not self.initiated:
            raise NotImplementedError(
                    'no init method called before launch')
        self.initFunctions()

        dt1 = dt_close *speed_factor[0]
        dt2 = dt_ip    *speed_factor[1]
        dt3 = dt_close *speed_factor[2]
        start = self.start
        stop  = self.start + tN
        max_steps   = (stop-start)/(dt2) + 1500001 
        time_passed = start
        dt  = dt1
        dtt = dt*dt

        pos0 = self.pos0
        vel0 = self.vel0
        self.start_planet = np.argmin( 
                np.linalg.norm(self.posFunc(start),axis=0)
                )
        self.target_planet = target_planet

        times = np.zeros((max_steps+1))
        sat_pos = np.zeros((max_steps+1,2))
        sat_vel = np.zeros_like(sat_pos)
        sat_acc = np.zeros_like(sat_pos)
        times[0] = time_passed
        sat_pos[0] = pos0
        sat_vel[0] = vel0
        sat_acc[0] = self.accelFunc(sat_pos[0], time_passed)
        i = 0
        
        prev_closest = 110 #random high value (AU)
        slow = True
        norm = np.linalg.norm
        break_stop = False
        
        while time_passed <= stop:
            if i >= max_steps or break_stop:
                if not break_stop:
                    print "Encountered max index! Breaking (perhaps) early"
                else:
                    print "Breaking manually"
                break
            sat_pos[i+1] = sat_pos[i] + sat_vel[i]*dt + 0.5*sat_acc[i]*dtt
            #staggered lf: temp_vel = sat_vel[i] +sat_acc[i-1]*0.5*dt
            sat_acc[i+1] = self.accelFunc(sat_pos[i+1], time_passed)
            sat_vel[i+1] = sat_vel[i] + 0.5*(sat_acc[i] + sat_acc[i+1])*dt
            #euler-cromer:
            #sat_vel[i+1] = sat_vel[i] + sat_acc[i]*dt
            #sat_pos[i+1] = sat_pos[i] + sat_vel[i+1]*dt
            #sat_acc[i+1] = self.accelFunc(sat_pos[i+1], time_passed)

            time_passed += dt
            times[i+1] = time_passed
            i+=1
            if i%1000 ==0:
                rel_dist = norm(self.posFunc(time_passed).T-sat_pos[i],axis=1)
                closest = np.amin(rel_dist)
                min_index = np.argmin(rel_dist)
                acc_now = norm(sat_acc[i])
                vel_now = norm(sat_vel[i])
                print "dist to cl(%d): %9f | a%8.2f | v%5.2f | t-left: %6.3f"\
                    %(min_index,closest,acc_now, vel_now, stop-time_passed)
                if min_index == target_planet and break_close:
                    if closest > prev_closest:
                        print 'found minimum', closest
                        break
                    prev_closest = closest
                force_relation = norm(sat_pos[i])*np.sqrt(self.starMass*self.mass)

                if slow:
                    '''Sets timestep larger/slower where force from 
                    closest planet is as large as force from star times
                    factor k = 2
                    '''
                    if np.all(rel_dist[min_index] \
                              > force_relation[min_index]/2.):
                        print "i = ", i
                        print "time: ", time_passed
                        dt = dt2
                        dtt = dt*dt
                        print "dt = ", dt
                        slow = False
                else:
                    if np.any(rel_dist[min_index] \
                              < force_relation[min_index]/2.): 
                        dt =  dt3
                        dtt = dt*dt
                        print "dt = ", dt
                        slow = True
                        

        # --- end of while loop ---

        self.sat_pos = sat_pos[:i]
        self.sat_vel = sat_vel[:i]
        self.sat_acc = sat_acc[:i]
        self.times = times[:i]

        self.initiated = False
        return sat_pos[:i], sat_vel[:i], sat_acc[:i], times[:i]

    def load_sim(self, fname, folder):
        "fname should be list of names of pos, vel and times"
        #filenames = 
        self.sat_pos = np.load( folder + fname[0] )
        self.sat_vel = np.load( folder + fname[1] )
        self.times = np.load( folder + fname[2] )

    def plot(self, planet = 5, end_at_closest_approach = False,
             plot_length= 5000):
        try:
            pos = self.sat_pos
            vel = self.sat_vel
            times = self.times
        except:
            print "simulation must be run before plotting"
            import sys
            sys.exit()

        inner_planets = np.array([0,1,4,5])
        planet_pos = self.posFunc(times)
        planet_vel = self.velFunc(times)

        relpos = pos - planet_pos[:,planet].T
        relvel = vel - planet_vel[:,planet].T
        rel_r = np.linalg.norm(relpos,axis = 1)
        rel_v = np.linalg.norm(relvel,axis = 1)
        closest_i = np.argmin(rel_r)
        print "Closest encounter with %d" %planet, rel_r[closest_i]
        print "Rel speed at CE with %d: " %planet, rel_v[closest_i]

        #relpos *= s.au /1000.
        plt.scatter(0,0,c='y')
        if  end_at_closest_approach:
            end = closest_i
        else:
            end = -1
        plt.plot(pos[:end,0], pos[:end,1])
        plt.plot(planet_pos[0,inner_planets, :end].T,
                planet_pos[1,inner_planets, :end].T)
        plt.legend(('satelite', 'planet0'))
        plt.axis('equal')
        plt.show()

        theta = np.linspace(0,2*np.pi,20)
        x = self.radius[planet] * np.cos(theta)
        y = self.radius[planet] * np.sin(theta)
        km_relpos = 1./1000 * self.au * relpos
        plt.plot(km_relpos[-plot_length:,0], km_relpos[-plot_length:,1])
        plt.scatter(km_relpos[-1,0],km_relpos[-1,1],s=50,c='g')
        plt.scatter(km_relpos[closest_i,0], km_relpos[closest_i,1], c='k')
        plt.plot(x,y)
        plt.axis('equal')
        plt.show()
        return rel_r[closest_i], closest_i

    def saveData(self, fname = ('satPos', 'satVel','satTimes'),
            folder='data'):
        try:
            pos = self.sat_pos
            vel = self.sat_vel
            times = self.times
            data = [pos,vel,times]
        except:
            print "simulation must be run before saving"
            import sys
            sys.exit()

        import glob
        num_old = len(glob.glob('%s/%s*'%(folder,fname[0])))
        for i,name in enumerate(fname):
            full_name = '%s/%s%d' %(folder,name,num_old)
            np.save(full_name, data[i])
            print "Saved file ", full_name


    def dragForce(self, pos, vel):
        R = self.radius_target_SI
        r = np.linalg.norm(pos)
        h =  r - R

        M = self.mass_target_SI
        G = self.G_SI
        rho = self.rhoFunc
        r = np.linalg.norm(pos)
        rel_vel = vel - self.atmosphere_vel(pos)
        v = np.linalg.norm(rel_vel)
        A = self.area
        C_D = 1
        return - 0.5*rho(h)*C_D*A*rel_vel*v

    def gravForce(self, pos):
        M = self.mass_target_SI
        m = self.sateliteMass
        G = self.G_SI
        r = np.linalg.norm(pos)
        return -pos*M*m*G/r**3

    def atmosphere_vel(self, pos):
        T = self.period[5]*self.day2sec #s
        circum = 2*np.pi*np.array((-pos[1],pos[0],0))
        return circum/T


    def orbit_sim(self, tN = 3600., dt = 20.,target_height=70e3):
        '''
        tN : float
            years to run sim
            
        dt : float
            seconds for each timestep
        '''
        
        from scipy.constants import G
        radius = self.radius[5] * 1000
        M = self.mass[5] * 1.98855e30

        sendSatelite = False
        if sendSatelite:
            t0= 21.3817592525
            position = (1.68932738, -5.98807047, 0)
            velocity = (3.16335623,  0.32657798, 0)
            star_pos0 = np.array(position)
            star_vel0 = np.array(velocity)
            posFunc = self.getPositionsFunction(xyz=True)
            velFunc = self.getVelocityFunction()
            planetPos = posFunc(t0)
            planetVel = velFunc(t0)
            pos0 = (planetPos[:,5] - star_pos0) * self.au 
            vel0 = (planetVel[:,5] - star_vel0) * self.au / self.year
        else:
            pos0 = np.array(( -6582268.62415 ,  -1141336.5617 ,  0 ))
            vel0 = np.array(( 439.800144275 ,  -2551.96645169 ,  0 ))

        target_orbit = target_height + radius #m
        rA = np.linalg.norm(pos0)
        rB = target_orbit
        at = (rA+rB)/2.  # semi major axis of transfer orbit
        vA = np.sqrt(G*M*(2./rA - 1./at)) # velocity wanted
        tTransfer = np.pi*np.sqrt((rA+rB)**3/(8*G*M))
        self.tOrbit = 2*tTransfer
        boost = - vel0 + np.array((0,0,vA))
        vel0 = vel0 + boost
        print boost

        vB = np.sqrt(G*M/rB)
        
        # --------------------------------------
        start = 0
        stop  = tN 
        #raw_input('continue?')
        n   = int((stop-start)/(dt)+1)

        self.radius_target_SI = radius
        self.mass_target_SI =  M
        self.rho0_target = self.rho0[5]
        self.G_SI = G
        self.rhoFunc = self.getDensityFunction()
        m = 1190  #kg
        self.sateliteMass = m
        self.area = 15 #6.2 #m^2
        dragForce = self.dragForce; gravForce = self.gravForce
        times = np.linspace(start,stop,n+1)
        pos = np.zeros((n+1,3))
        vel = np.zeros_like(pos)
        F_d = np.zeros_like(pos)
        F_g = np.zeros_like(pos)

        times[0] = start
        pos[0] = pos0
        vel[0] = vel0
        F_d[0] = dragForce(pos[0],vel[0])
        F_g[0] = gravForce(pos[0])
        acc =  (F_d[0] + F_g[0] )/m

        print acc
        for i in range(n):
            vel[i+1] = vel[i] + acc*dt
            pos[i+1] = pos[i] + vel[i+1]*dt
            F_d[i+1] = dragForce(pos[i+1],vel[i+1])
            F_g[i+1] = gravForce(pos[i+1])
            acc =  (F_d[i+1] + F_g[i+1] )/m
            if np.linalg.norm(pos[i+1]) < radius:
                print "Inside planet"
                break
            if i%((n)/100) == 0:

                print "%.1f%%" % (float(1+100*i)/n)
                #print acc
                #print F_g[i+1]
                #print np.linalg.norm(F_d[i+1])
                #raw_input()
                
        self.pos = pos[:i]
        self.vel = vel[:i]
        self.F_d = F_d[:i]
        self.F_g = F_g[:i]
        self.acc = (F_d[:i] + F_g[:i] )/m
        self.times = times[:i]
        return self.pos, self.vel, self.acc, self.times

        

# 7.45955363121e-05 collision t14.25, a1.67pi, v9590

if __name__ == "__main__":
    s = MySateliteSim(87464)
    s.launch_init()
    s.satelite_sim()
    s.plot()
    if raw_input('save?(y/N)') == 'y':
        b = 9594
        fnames = [w+b for w in ['pos','vel','times']]
        s.savedata(fname = fnames)

