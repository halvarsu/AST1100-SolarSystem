from MySateliteSim import MySateliteSim
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def rot_matrix(angle):
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array(((c,-s),(s,c)))

def boost(velocity, boost, boost_angle):
    direction = velocity/np.linalg.norm(velocity)
    boost_direction = rot_matrix(boost_angle).dot(direction)

    new_vel = vCl + boost_direction*boost/system.au*system.year
    return new_vel



system = MySateliteSim(87464)
au = system.au
launch_v = 9594

p = np.load('positionsHomePlanet.npy')
t = np.load('times.npy')
st = np.load('data/times%d.npy' %launch_v)
s = np.load('data/pos%d.npy' % launch_v)

posFunc = interp1d(t[:-1], p)
planetPos = posFunc(st)

satFunc = interp1d(st, s.T)
def vel(t,dt=0.0000001):
    return (satFunc(t+dt)-satFunc(t-dt))/(2*dt)

planet4 = planetPos[:,4]
relPos = s.T-planet4
relDist = np.linalg.norm(relPos, axis =0)
closest_dist = np.amin(relDist[:-10])
closest_i = np.argmin(relDist[:-10])
relposCl = relPos[:,closest_i]
hCl = closest_dist - system.radius[4]*1000./au
tCl = st[closest_i]
pCl = s[closest_i]
aCl = np.arctan(relposCl[1]/relposCl[0])#angle
vCl = vel(tCl)

#if hCl < 0:
    #print "Wont bother with values inside planet"
    #continue

print "time of closest approach: ", tCl
print "satpos at closest approach: ", pCl
print "distance from surface at closest approach: ", hCl
print "relpos at closest approach", relposCl
print "angle at closest approach", aCl

boost_v = 2000

boost_a = 1.07*np.pi/4
vel0 = boost(vCl, boost_v, boost_a)
pos0 = pCl 
sim_length = 6.2
close_speed_factor = 1

# 0.336138702386  b1000 a1.06pi/4 lv 9594
# 0.0596081065397 b2000
# 0.00503893068436 a1.07pi/4
# 0.0268986034289 a1.065pi/4 v0.629758507124 


dt_far = 1/(50000.)
print "---------------------------------------------"
print "continuing the launch ", launch_v
psim, vsim, asim, tsim = system.satelite_sim(t0 = tCl, 
        tN = tCl+sim_length, start_planet=4, target_planet =5,
        vel0 = vel0, pos0 = pos0, first_boost = False,
        dt_far=dt_far, sec_per_dt_close = close_speed_factor)
filenames = ['%sv%db%d' %(name,launch_v,boost_v) for name in ('pos','vel','times')]
system.saveData(fname = filenames, folder = 'postGAdata')
system.plot()

