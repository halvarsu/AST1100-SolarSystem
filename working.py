# coding: utf-8
from MySolarSystem import MySolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def simulate(seed, System, years = 8,planet_index = 0, units='au'
        , res = 10000):
    au =  149597870700          #m
    sec_in_year = 365 * 24 * 3600      #s
    G = 6.67408e-11 
    solar_mass = 1.98855e30     #kg
    n = res *years
    stop = years * sec_in_year 
    dt = stop/n+1

    starMass_kg = System.starMass * solar_mass
    mass_kg = System.mass[planet_index] * solar_mass

    x0 = System.x0 * au
    y0 = System.y0 * au
    vx0 = System.vx0 * au / sec_in_year
    vy0 = System.vy0 * au / sec_in_year

    times = np.linspace(0,stop,n)
    pos = np.zeros((n+1, 2))
    vel = np.zeros((n+1, 2))

    pos[0] = np.array((x0[planet_index], y0[planet_index]))
    vel[0] = np.array((vx0[planet_index], vy0[planet_index]))
    accel = np.zeros((n+1, 2))
    #accel[0] = calculate_forces(pos[0]) / mass_kg # second planet

    r = np.sqrt(pos[0][0]**2 + pos[0][1]**2)
    accel[0] = -G*starMass_kg*pos[0]/r**3

    for i in range(len(times)): 
        r = np.sqrt(pos[i][0]**2 + pos[i][1]**2)
        accel[i] = -G*starMass_kg*pos[i]/r**3
        #accel[i] =  calculate_forces(pos[i]) / mass_kg

        vel[i+1] = vel[i] + accel[i] * dt
        pos[i+1] = pos[i] + vel[i+1] * dt
        #print accel[i], vel[i]

    if units == 'au':
        return pos/au
    else:
        return pos

def sim2(seed, System, years = 8, planet_index = "", res = 10000):
    G = 4*np.pi**2
    n = res*years
    stop = years
    dt = stop/n+1

    starMass = System.starMass
    mass = System.mass





def planet_plotter(positions, years, fig=None, ax = None):
    if not ax:
        fig, ax = plt.subplots()
    n = len(positions)
    ax.plot(positions.T[0], positions.T[1])
    colors = np.linspace(0,1,years+1)
    ax.scatter(positions.T[0][::n/years], positions.T[1][::n/years], c=colors)
    ax.axis('equal')
    return fig, ax
    

if __name__ == "__main__":
    seed = 87464
    MyStarSystem = MySolarSystem(seed) 
    years = 10
    fig, ax = plt.subplots()
    ax.scatter(0,0, c='y')

    for x in range(MyStarSystem.numberOfPlanets):
        print x
        positions = simulate(seed, MyStarSystem, planet_index = x
                ,years=years, res = 5000)
        fig, ax = planet_plotter(positions, years, fig, ax)

    plt.axis('equal')
    plot_max = max(MyStarSystem.a +1)
    plt.axis([-plot_max,plot_max,-plot_max,plot_max])
    ax = MyStarSystem.polar_orbits(ax =ax)
    plt.show()
