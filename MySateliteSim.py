# coding: utf-8
from MySolarSystem import MySolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings


ang2pix = MySolarSystem.ang2pix
import warnings
warnings.filterwarnings("ignore")


def rot_matrix(angle):
    """Provides a two dimensional rotation matrix for a given angle"""
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array((c,-s),(s,c))



class MySateliteSim(MySolarSystem):
    """This is an implementation of My Solar System designed for sending
    satelites into space"""
    au = 149597870700 # m
    year = int(365.25*24*3600) #s

    def __init__(self, seed): 
        """
        On init, creates an instance of MySolarSystem with the seed
        provided, and loads arrays with planet positions and times. 
        These are then used to create interpolation functions for 
        position, velocity and acceleration of the planets.

        :seed: The seed to be instantiated in MySolarSystem

        """
        MySolarSystem.__init__(self, seed)

        self._seed = seed
        self.planetPos  =  np.load('positionsHomePlanet.npy')
        self.times      =  np.load('times.npy')
        self.posFunc    =  self.getPositionsFunction()
        self.velFunc    =  self.getVelocityFunction()
        self.angleFunc  =  self.getAngleFunction()
        self.accelFunc  =  self.getAccelFunction()


    def testVelFunction(self):
        velFunc= self.getVelocityFunction()
        exact = np.linalg.norm((self.vx0[0], self.vy0[0]))
        success = True
        for i in range(10):
            test = np.linalg.norm(velFunc(0, 10**-i), axis = 0)[0]
            print test, abs((test-exact)/exact), "   ", 10**-i
        print exact

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
        if not hasattr(self, 'accelFunction'):
            self.accelFunction = np.vectorize(self.accelFunc_unvec)
        return self.accelFunction

    def boost(self, vel, boost, boost_angle):
        """
        Applies an instantanious boost at an angle boost_angle at the 
        vel direction

        :vel:      two dimensional velocity before boost
        :boost:         magnitude of boost
        :boost_angle:   angle of boost
        """
        direction = vel/np.linalg.norm(vel)
        boost_direction = np.dot(rot_matrix(boost_angle), direction)
        new_vel = vel + boost_direction*boost/self.au*self.year
        return new_vel

    def boost_init(
            self, pos0, vel0, t0, init_boost = 0, boost_angle = 0,
            target_planet = 4):
        self.start = t0
        self.vel0 = self.boost(vel0, init_boost, boost_angle)
        self.pos0 = pos0
        self.start_planet = np.argmin(
                                np.linalg.norm(self.posFunc(t0),axis=0))
        self.target_planet = target_planet
        return self.start_planet, self.vel0

    def launch_init(
            self, t0 = 0, init_boost = 10000, init_angle = 1.6*np.pi, 
            start_dist=1000, start_planet = 0, target_planet = 4,
            ground_launch = True):
        self.start = t0
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
            init_pos =  planet_direction * self.radius[start_planet]
            init_vel = planet_direction *init_boost
        else:
            # init_angle important for direction of hyperbolic excess vel
            theta = init_angle + planet_angle 
            init_vel = np.array((-init_boost*np.sin(theta),
                                    init_boost*np.cos(theta)))
            init_dist = start_dist+self.radius[start_planet]
            init_pos = np.array((init_dist*np.cos(theta),
                                 init_dist*np.sin(theta)))

        pos0 = planetPos0.T[start_planet] + init_pos/self.au*1000
        vel0 = planetVel0.T[start_planet] + init_vel*self.year/self.au
        self.pos0 = pos0
        self.vel0 = vel0
        self.start_planet  = start_planet 
        self.target_planet = target_planet
        return self.start_planet, self.vel0


    def satelite_sim(
        self, tN = 4, dt_ip = 1/(50000.), dt_close = 1/(365.25*24*3600), 
        speed_factor = [1,1,1], *args, **kwargs):

        ''' 
        Simulates satelite trajectory . There are two options, either a
        first_launch from planet, or a continuation from a point in space.

        First 
        continuing from a point in space if first_launch = False.
        ground_launch
        The satelite is launched using 
        Accuracy goes down with higher speed_factor
        [tN] = years
        ''' 
        
        dt1 = dt_close *speed_factor[0]
        dt2 = dt_ip    *speed_factor[1]
        dt3 = dt_close *speed_factor[2]
        start = self.start
        stop  = tN
        max_steps   = (stop-start)/(dt2) + 1500001 
        
        time_passed = start
        dt  = dt1
        dtt = dt*dt

        pos0 = self.pos0
        vel0 = self.vel0
        start_planet = self.start_planet
        target_planet = self.target_planet


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
                print "relative dist to closest(%d): %f  acc: %f  t left: %.3f"\
                        %(min_index,closest,acc_now, tN-time_passed)

                if min_index == target_planet:
                    if closest > prev_closest:
                        print 'found minimum', closest
                        break
                    old_closest = closest
                force_relation = norm(sat_pos[i])*np.sqrt(self.starMass*self.mass)

                if slow:
                    '''Sets timestep larger/slower where force from closest
                    planet is as large as force from star
                    '''
                    if np.all(rel_dist[min_index]>force_relation[min_index]/2.):
                        print "i = ", i
                        print "time: ", time_passed
                        dt = dt2
                        dtt = dt*dt
                        print "dt = ", dt
                        slow = False
                        
                else:
                    if np.any(rel_dist[min_index]<force_relation[min_index]/2.): 
                        dt =  dt3
                        dtt = dt*dt
                        print "dt = ", dt
                        slow = True
                        

        # --- end of while loop ---

        self.sat_pos = sat_pos[:i]
        self.sat_vel = sat_vel[:i]
        self.sat_acc = sat_acc[:i]
        self.times = times[:i]

        return sat_pos[:i], sat_vel[:i], sat_acc[:i], times[:i]

    def load_sim(self, fname, folder):
        "fname should be list of pos, vel and times"
        #filenames = 
        self.sat_pos = np.load( folder + fname[0] )
        self.sat_vel = np.load( folder + fname[1] )
        self.times = np.load( folder + fname[2] )

    def plot(self, planet = 5):
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
        plt.plot(pos[:,0], pos[:,1])
        plt.plot(planet_pos[0,inner_planets].T,
                planet_pos[1,inner_planets].T)
        plt.legend(('satelite', 'planet0'))
        plt.axis('equal')
        plt.show()

        
        theta = np.linspace(0,2*np.pi,20)
        x = self.radius[planet] * np.cos(theta)
        y = self.radius[planet] * np.sin(theta)
        km_relpos = 1./1000 * self.au * relpos
        plt.plot(km_relpos[-5000:,0], km_relpos[-5000:,1])
        plt.scatter(km_relpos[-1,0],km_relpos[-1,1])
        plt.scatter(km_relpos[closest_i,0], km_relpos[closest_i,1])
        plt.plot(x,y)
        plt.axis('equal')
        plt.show()
        return rel_r[closest_i]
        
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


    def orient(self, time):
        planet_pos = self.posFunc(time)
        #rel_pos = 
        
        #for x,y in 
        

# 7.45955363121e-05 collision t14.25, a1.67pi, v9590

if __name__ == "__main__":
    s = MySateliteSim(87464)
    start_time = 14.25 #years
    stop_time = start_time + 3. #years
    init_boost = 9587#m/s
    angle = 1.67*np.pi #radians
    import time
    sim_start = time.time()
    dt_close = 1/(365.25*24*3600)
    for x in range(3):
        s.launch_init(t0 = start_time, 
                init_boost = init_boost-x,init_angle = angle,
                ground_launch = True)
        pos, vel, acc, times= s.satelite_sim(
            t0=start_time, init_boost=init_boost-x, tN = stop_time, 
            init_angle = angle, dt_close = dt_close
        )
        print "length of pos ", pos.shape[0]
        print "Time used on sim: ", time.time() - sim_start,"s"
        norm = np.linalg.norm
        print "Last value: ", norm(pos[-1])
        print "Distance covered: ", norm(pos[0]-pos[-1])
        #s.plot()
        filenames = ['%s%d' %(name,init_vel-x) for name in ('vel','pos','times')]
        s.saveData(filenames)
