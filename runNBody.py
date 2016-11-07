import matplotlib.pyplot as plt
import numpy as np
import sys
from NBodySystem import NBodySystem
import argparse
import seaborn as sns

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--sim', action = 'store_true', default=False,
                    help='Simulate instead of loading existing files')
parser.add_argument('-y','--years', type=int,default=200)
parser.add_argument('-r','--plot_res',
                    help='Resolution of plotting',type=int, default=1)
parser.add_argument('-p','--plot_nr',
                    help='Which of the plots to produce',type=int,
                    default=1,choices=range(1,6))
args = parser.parse_args()

seed = 87464
system = NBodySystem(seed)
max_years = 200; res = 2000
years = args.years

if args.sim:
    pos, vel, acc, times = system.simulate_nbody(years = years, res = res)
else:
    pos, vel, acc, times = system.loadData()
    if years < max_years:
        pos = pos[:years*res]; vel = vel[:years*res]
        acc = acc[:years*res]; times = times[:years*res]
    elif years > max_years:
        print "Warning, too many years. set to max = %d" %max_years
        years = max_years

pos = pos[::args.plot_res]
vel = vel[::args.plot_res]
acc = acc[::args.plot_res]
times = times[::args.plot_res]


#
if args.plot_nr == 1:
    plt.axis('equal')
    plot_max = max(system.a +1)
    plt.axis([-plot_max,plot_max,-plot_max,plot_max])
    plt.plot(pos[:,0], pos[:,1])
    plt.show()

if args.plot_nr == 2:
    plt.axis('equal')
    plt.plot(pos[:,0,-1], pos[:,1,-1])
    plt.title('Star movement over %d years' %years)
    plt.savefig('figure/starMovement.png')
    plt.show()

if args.plot_nr == 3:
    fig, axes = plt.subplots(3, sharex=True)
    axes[0].plot(times, pos[:,:,-1])
    axes[0].set_title('Star position')
    axes[1].plot(times, vel[:,:,-1])
    axes[1].set_title('Star velocity')
    axes[2].plot(times, acc[:,:,-1])
    axes[2].set_title('Star acceleration')
    prefix = ['','v_','a_']
    for pre, ax in zip(prefix, axes):
        legend = ax.legend(["$%sx$"%pre,"$%sy$"%pre],fontsize=12)
    plt.savefig('figure/nBodyStar.png')
    plt.show()

if args.plot_nr == 4:
    ax = plt.subplot(111)
    v_star = np.zeros((len(vel),3))
    v_star[:,:2] = vel[:,:,-1]
    i = np.pi/4
    v_pec = 2 #AU/yr
    los = np.pi/4
    uvec = np.array((np.cos(i)*np.cos(los),np.cos(i)*np.sin(los),np.sin(i)))
    v_rstar = np.einsum('ij,j->i',v_star,uvec)


    mu = 0; sigma = 0.2*np.max(np.abs(vel[:,0,-1]))
    noise = np.random.normal(mu,sigma,len(times))

    ax.plot(times, v_rstar +v_pec+ noise,linewidth=0.07)
    ax.set_title('$i=\pi/4$, $v_{pec} = 2 \\frac{AU}{yr}$,'\
            ' $\\theta=\pi/4$ ',fontsize = 14)

    ax.set_xlabel('years', fontsize=12)
    ax.set_ylabel('Peculiar velocity $v_{r*}$', fontsize = 12)

    plt.savefig('figure/starVelWithNoise.png')
    plt.show()

if args.plot_nr == 5:
    starPos = pos[:,:,-1]
    planetPos = pos[:,:,:-1]
    
    flux = np.ones(len(times),dtype='float')
    mu = 0; sigma = 0.2
    noise = np.random.normal(mu,sigma,len(times))

    radius = system.radius *1000/ float(system.au)
    starRadius = system.starRadius *1000/float(system.au)
    mask = planetPos[planetPos!=planetPos]

    for i in range(system.numberOfPlanets):
        cover = radius[i]**2/starRadius**2
        pPos = planetPos[:,:,i]
        mask1 = np.abs(pPos[:,1]-starPos[:,1])<(starRadius+radius[i])
        print cover
        mask2 = pPos[:,0] > 0
        mask3 = np.logical_and(mask1,mask2)
        flux[mask3] -= cover
        if i == 2:
            index = np.argsort()[-1]
            

    print radius**2/starRadius**2
    plt.plot(times*365,flux+noise)
    plt.scatter(times[index], 0.8)
    plt.show()

    print


    # line of sight = x-axis


