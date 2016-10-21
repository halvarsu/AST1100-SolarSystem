import PyKEP as pkp
from PyKEP import epoch, AU, MU_SUN, DEG2RAD, lambert_problem, DAY2SEC, DAY2YEAR
from PyKEP.orbit_plots import plot_planet, plot_lambert
import numpy as np
from MySateliteSim import MySateliteSim
import matplotlib as mpl
import matplotlib.pyplot as plt




system = MySateliteSim(87464)
def init_planet(system, planet = 0, name = 'none'):
    if name == 'none':
        name = str(planet)
    SEC2YEAR = DAY2SEC/float(DAY2YEAR)
    x0 = system.x0[planet]*AU
    y0 = system.y0[planet]*AU
    vx0 = system.vx0[planet]*AU/SEC2YEAR
    vy0 = system.vy0[planet]*AU/SEC2YEAR
    mu_star = system.starMass * MU_SUN
    radius = system.radius[planet] * 1000
    mu_planet = system.mass[planet] * MU_SUN
    safe_radius = radius * 1.2
    planet = pkp.planet.keplerian(
            epoch(0), np.array((x0,y0,0)), np.array((vx0,vy0,0)),
            mu_star, mu_planet, radius, safe_radius, name)
    return planet


def init_planet_old2(system, planet=0, name = 'home'):
    a,e,i,W,w,M = system.getKeplerianElements(new_values = True).T[planet]
    mu_star = system.starMass * MU_SUN
    mu_planet = system.mass[planet] * MU_SUN 
    radius = system.radius[planet] * 1000
    safe_radius = radius + 100000
    return pkp.planet.keplerian(epoch(0), np.array((a*AU, e, i,
            W, w, M)), mu_star, mu_planet, radius, safe_radius,name)

def init_planet_old(planet):
    x0 = system.x0[planet]
    y0 = system.y0[planet] 
    e = system.e[planet]
    a = system.a[planet]
    b = a*np.sqrt(1-e**2)
    M = np.arctan((y0*a)/(x0*b)) - e*y0/b + np.pi
    i = 0.0
    W = 0.0
    w = system.psi[planet]  + np.pi
    mu_star = system.starMass * MU_SUN
    mu_planet = system.mass[planet] * MU_SUN 
    radius = system.radius[planet] * 1000
    safe_radius = radius + 100000
    return pkp.planet.keplerian(epoch(0), np.array((a*AU, e, i,
            W, w, M)), mu_star, mu_planet, radius, safe_radius, 'home')

def lambert_solve(t1,t2,A,B, nameA='herbin', nameB='hars'):
    t1 = epoch(t1/DAY2YEAR)
    t2 = epoch(t2/DAY2YEAR)
    dt = (t2.mjd2000 - t1.mjd2000)*DAY2SEC

    planet1 = init_planet(system,A, name = nameA)
    planet2 = init_planet(system,B, name = nameB)

    r1, v1 = planet1.eph(t1)
    r2, v2 = planet2.eph(t2)
    #print np.array(r1) / AU

    mu_star = system.starMass * MU_SUN
    l = lambert_problem(r1, r2, dt, mu_star)
    return l


def lambert_plot(t1,t2, A, B):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    axis = fig.gca(projection = "3d")

    t1 = epoch(t1/DAY2YEAR)
    t2 = epoch(t2/DAY2YEAR)
    dt = (t2.mjd2000 - t1.mjd2000)*DAY2SEC

    planet1 = init_planet(system, A)
    plot_planet(planet1, t0=t1, legend=True, units = AU, ax = axis)
    planet4 = init_planet(system, B)
    plot_planet(planet4, t0=t2, legend=True, units = AU, ax = axis)

    rH, vH = planet1.eph(t1)
    rT, vT = planet4.eph(t2)

    mu_star = system.starMass * MU_SUN
    l = lambert_problem(rH, rT, dt, mu_star)
    plot_lambert(l, legend=True, units = AU, ax = axis)
    #print j
    #plot_lambert(l, sol=1, legend=True, units = AU, ax = axis)
    axis.set_xlabel('x')
    plt.show()
    return l

def test_init_and_solve():
    system = MySateliteSim(87464)
    year = system.year
    for x in range(1,3):
        planet = init_planet(system, planet = x)

        eph = np.array(planet.eph(epoch(0)))
        print eph/AU
        print "-------------------------"
        print  x
        print "Position:"
        print "from pykep  x %7.4f  y %7.4f  " %(eph[0,0]/AU, eph[0,1]/AU)
        print "from module x %7.4f  y %7.4f  " %(system.x0[x], system.y0[x])
        print "Velocity:"
        print "from pykep  vx %7.4f  vy %7.4f  " %(eph[1,0]/AU*year,
                                                   eph[1,1]/AU*year)
        print "from module vx %7.4f  vy %7.4f  " %(system.vx0[x], system.vy0[x])

    t1 =16.9105290335 
    t2 =20.6651033379
    A = 4
    B = 5

    l = lambert_solve(t1,t2,A,B)
    lambert_plot(t1,t2,A,B)

def plot_all(t = 0):
    from mpl_toolkits.mplot3d import Axes3D
    system = MySateliteSim(87464)
    fig = plt.figure()
    axis = fig.gca(projection = "3d")
    t0 = epoch(t/DAY2YEAR)
    for x in range(system.numberOfPlanets):
        Planet = init_planet(system, planet = x)
        plot_planet(Planet, t0=t0, legend=True, units = AU, ax = axis)
    plt.show()





if __name__ == "__main__":
    #lambert_solve(14.32, 17.32, 0 ,4)
    plot_all()
