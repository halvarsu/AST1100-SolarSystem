import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter

def find_angle(pos1, pos2):
    '''Expects two tuples or arrays with 2 as first shape'''
    a1 = np.arctan(pos1[1], pos1[0])
    a2 = np.arctan(pos2[1], pos2[0])
    return a1 - a2


def find_hohmann_params(system, A=0, B=1, t0 = 0):
    '''returns parameters of the hohmann transfer from planet A to planet B
    at time t0 after start of sim (epoch)
    
    :system: : Expects instance of class MySolarSystem
    '''
    G = 4*np.pi**2
    i1 = A
    i2 = B

    planetPos = np.load('positionsHomePlanet.npy')
    times = np.load('times.npy')

    posFunction = inter.interp1d(times[:-1], planetPos)
    theta, orbits_r = system.analytical_pos(planets = "%s, %s"% (A,B))
    print theta.shape, orbits_r.shape
    angleFunction = inter.interp1d(theta, orbits_r.T)

    # Analytical_pos returns positions as function of theta, while
    # posFunction returns positions as function of time.

    orbit_B_r = orbits_r.T[B]

    posA = posFunction(t0).T[A]
    rA = np.linalg.norm(posA)
    print rA 
    angleB = find_angle(posA, (0,0)) + np.pi # The angle of pos of planet B
    t_index = int(len(theta)*angleB)         # index used to find pos
    rB = orbit_B_r[t_index]
    #rB = np.norm(posFunction(t0+system.period[B]/2)[B])

    aTransfer = (rA+rB)/2.
    GM = G * system.starMass
    v_init_A     = np.sqrt(GM/rA)
    v_final_B    = np.sqrt(GM/rB)
    v_transfer_A = np.sqrt(GM*(2/rA - 1/aTransfer))
    v_transfer_B = np.sqrt(GM*(2/rB - 1/aTransfer))
    delta_VA = v_init_A - v_transfer_A
    delta_VB = v_final_B - v_transfer_B
    total_delta_V = delta_VA + delta_VB
    e = 1 - rA/aTransfer #eccentrisity of transfer ellipse
    transferTime = np.pi*np.sqrt((rA+rB)**3/(8*GM))
    return transferTime
    
if __name__ == "__main__":
    import MySolarSystem as mySS
    system = mySS.MySolarSystem(87464)

    print find_hohmann_params(system)
