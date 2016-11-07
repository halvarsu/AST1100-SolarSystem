from MySateliteSim import MySateliteSim 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
system = MySateliteSim(87464)

parser = argparse.ArgumentParser()
parser.add_argument('-p','--part', help='choose plot-part. -1>all,0>none',
                    default = -1, choices=range(5),type=int)
parser.add_argument('-l','--load', action='store_true', default=False ,
                    help='load values instead of simulating')


args = parser.parse_args()

sim_days = 0.6
if not args.load:
    p,v,a,t = system.orbit_sim(tN = 3600*24*sim_days, dt = 0.5, 
            orbit_height= 320e4)
    np.save('data/tOrbit', system.tOrbit)
    np.save('data/orbit_pos', p)
    np.save('data/orbit_vel', v)
    np.save('data/orbit_acc', a)
    np.save('data/orbit_times', t)
    F_d = system.F_d; F_g = system.F_g
    forces = np.array((F_d, F_g))
    np.save('data/orbit_forces', forces)
    np.save('data/orbit_boosts', (system.circ_i, system.land_i))
else:
    p = np.load('data/orbit_pos.npy')
    v = np.load('data/orbit_vel.npy')
    a = np.load('data/orbit_acc.npy')
    t = np.load('data/orbit_times.npy')
    F_d, F_g = np.load('data/orbit_forces.npy')
    system.circ_i, system.land_i = np.load('data/orbit_boosts.npy')


r = system.radius[5] * 1000
theta = np.linspace(0,2*np.pi, 100)
x = r * np.cos(theta)
y = r * np.sin(theta)

if args.part == 1 or args.part == -1:
    if 1:
        from plot_parametric import plot_parametric, sphere
        planet_sphere = lambda u,v : sphere(u,v,r=r)
        ax = plot_parametric(planet_sphere, return_ax = True)
    else:
        from plot_implicit import plot_implicit, sphere
        planet_sphere = lambda x,y,z: sphere(x,y,z,r=r)
        ax = plot_implicit(planet_sphere, bbox=(-r,r))

    plt.plot(p[:,0], p[:,1], zs = p[:,2])
    #plt.scatter(p[0,0],p[0,1],p[0,2])
    plt.show()
    #plt.plot(v[:,0], v[:,1])
    #plt.axis('equal')
    #plt.show()

if args.part == 2 or args.part == -1:
    plt.plot(x,y)
    #u_vec = v[0,:2]/np.linalg.norm(v[0,:2])

    i1 = system.circ_i
    i2 = system.land_i
    xy = np.linalg.norm(p[:,:2],axis=1)*np.sign(-p[:,0])
    z = p[:,2]
    plt.plot(xy, z)
    print ' circ_pos : ', p[i1]/np.linalg.norm(p[i1])
    print ' land_burn_pos : ', p[i2]/np.linalg.norm(p[i2])
    print ' land_pos : ', p[-1]/np.linalg.norm(p[-1])
    print ' lander init_vel :', v[i2]
    print ' lander launch time :', t[i2]
    m_pos = [-3341783.05623727, -598581.99217862,  289480.1740502 ]


    m_xy = np.linalg.norm(m_pos[:2])
    m_z = m_pos[2]
    plt.scatter(xy[i1], z[i1])
    plt.scatter(xy[i2], z[i2])
    plt.scatter(m_xy, m_z, c = 'g', s = 50)
    plt.axis('equal')
    plt.show()

if args.part == 3 or args.part == -1:
    vel = np.linalg.norm(v, axis = 1)
    pos = np.linalg.norm(p, axis = 1)
    plt.plot(t,pos)
    plt.plot((t[0],t[-1]),(r,r))
    plt.title('position')
    plt.show()
    plt.plot(t,vel)
    plt.title('vel')
    plt.show()

def derivate(t, func):
    dt = t[1]-t[0]
    dfdt = (func[1:] - func[0:-1])/dt
    return dfdt


if args.part == 4 or args.part == -1:
    fig, (ax1, ax2) = plt.subplots(2, sharex = True)
    F_d = np.linalg.norm(F_d, axis = 1)
    F_g = np.linalg.norm(F_g, axis = 1)
    ax1.plot(t,F_d)
    ax1.plot(t, F_g)
    ax1.legend(['drag', 'gravity'])
    ax2.plot(t[0:-1], derivate(t, F_d))
    ax3 = plt.subplot(111,projection='3d')
    h = (np.linalg.norm(p,axis = 1) - r)[-1000:]
    v = np.linalg.norm(v,axis = 1)[-1000:]
    h,v = np.meshgrid(h,v)
    plt.pcolor(h,v,F_d)
    ax3.plot(h,F_d)
    plt.show()



