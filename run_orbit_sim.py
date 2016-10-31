from MySateliteSim import MySateliteSim 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
system = MySateliteSim(87464)

parser = argparse.ArgumentParser()
parser.add_argument('-p','--part', default = 1, choices=range(5),type=int)
parser.add_argument('-l','--load', action='store_true', default=False )
args = parser.parse_args()

sim_days = 0.5
if not args.load:
    p,v,a,t = system.orbit_sim(tN = 3600*24*sim_days, dt = 1, 
            target_height= 128e3)
    np.save('data/tOrbit', system.tOrbit)
    np.save('data/orbit_pos', p)
    np.save('data/orbit_vel', v)
    np.save('data/orbit_acc', a)
    np.save('data/orbit_times', t)
else:
    p = np.load('data/orbit_pos.npy')
    v = np.load('data/orbit_vel.npy')
    a = np.load('data/orbit_acc.npy')
    t = np.load('data/orbit_times.npy')


r = system.radius[5] * 1000
theta = np.linspace(0,2*np.pi, 100)
x = r * np.cos(theta)
y = r * np.sin(theta)

if args.part == 1 or args.part == 0:
    #from plot_parametric import plot_parametric, sphere
    from plot_implicit import plot_implicit, sphere

    #planet_sphere = lambda u,v : sphere(u,v,r=r)
    #ax = plot_parametric(planet_sphere, return_ax = True, ax = ax)
    planet_sphere = lambda x,y,z: sphere(x,y,z,r=r)
    ax = plot_implicit(planet_sphere, bbox=(-r,r))
    plt.plot(p[:,0], p[:,1], zs = p[:,2])
    #plt.scatter(p[0,0],p[0,1],p[0,2])
    plt.show()
    #plt.plot(v[:,0], v[:,1])
    #plt.axis('equal')
    #plt.show()

if args.part == 2 or args.part == 0:
    plt.plot(x,y)
    #u_vec = v[0,:2]/np.linalg.norm(v[0,:2])

    plt.plot(np.linalg.norm(p[:,:2],axis=1)*np.sign(p[:,0]), p[:,2])
    plt.axis('equal')
    plt.show()

if args.part == 3 or args.part == 0:
    vel = np.linalg.norm(v, axis = 1)
    pos = np.linalg.norm(p, axis = 1)
    plt.plot(t,pos)
    plt.plot((t[0],t[-1]),(r,r))
    plt.title('position')
    plt.show()
    plt.plot(t,vel)
    plt.title('vel')
    plt.show()

if (args.part == 4 or args.part == 0 ) and not args.load:
    F_d = np.linalg.norm(system.F_d, axis = 1)
    F_g = np.linalg.norm(system.F_g, axis = 1)
    plt.plot(t,F_d)
    plt.plot(t, F_g)
    plt.legend(['drag', 'gravity'])
    plt.show()
