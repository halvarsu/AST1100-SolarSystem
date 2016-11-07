import numpy as np
import time
from MySolarSystem import MySolarSystem
from scipy.constants import m_p, m_e, k, G
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--only_isotermic', action = 'store_true',
                    default = False)

args = parser.parse_args()


if not args.only_isotermic:
    system = MySolarSystem(87464)
    system.pressure_solver()
else:
    sunMass = 1.989e30
    mu = 38
    m_H = m_p + m_e
    r_p = system.radius[5] * 1000
    M = system.mass[5] * sunMass
    rho0 = system.rho0[5] 
    T = 169

    def rho(r):
        return rho0 * np.exp((mu*m_H*G*M/(k*T))*(1/r - 1/r_p))

    r = np.linspace(r_p, r_p + 100000, 1000)

    rho_arr = rho(r)
    cent_r = np.argmin(np.abs((rho_arr - rho0/100)))
    print r[cent_r], r_p, r[cent_r] - r_p
    plt.plot(r-r_p, rho_arr)
    plt.scatter(r[cent_r]-r_p, rho_arr[cent_r])
    plt.show()
