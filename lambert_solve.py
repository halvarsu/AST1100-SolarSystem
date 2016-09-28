import PyKEP as pkp
from PyKEP import epoch, AU, MU_SUN, DEG2RAD, lambert_problem, DAY2SEC, DAY2YEAR
from PyKEP.orbit_plots import plot_planet, plot_lambert
import numpy as np
from MySateliteSim import MySateliteSim
import matplotlib as mpl
import matplotlib.pyplot as plt




system = MySateliteSim(87464)
def init_planet(planet):
    x0 = system.x0[planet]
    y0 = system.y0[planet] 
    e = system.e[planet]
    a = system.a[planet]
    b = a*np.sqrt(1-e**2)
    M = np.arctan((y0*a)/(x0*b)) - e*y0/b
    i = 0.0
    W = 0.0
    w = system.psi[planet]
    mu_star = system.starMass * MU_SUN
    mu_planet = system.mass[planet] * MU_SUN 
    radius = system.radius[planet] * 1000
    safe_radius = radius + 100000
    return pkp.planet.keplerian(epoch(0), np.array((a*AU, e, i,
            W, w, M)), mu_star, mu_planet, radius, safe_radius, 'home')

def lambert_solve(t1,t2,A,B):
    t1 = epoch(t1/DAY2YEAR)
    t2 = epoch(t2/DAY2YEAR)
    dt = (t2.mjd2000 - t1.mjd2000)*DAY2SEC

    planet1 = init_planet(A)
    planet2 = init_planet(B)

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

    planet1 = init_planet(A)
    plot_planet(planet1, t0=t1, legend=True, units = AU, ax = axis)
    planet4 = init_planet(B)
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

if __name__ == "__main__":
    lambert_solve(14.32, 17.32, 0 ,4)
