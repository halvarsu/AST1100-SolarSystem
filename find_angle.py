import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter

def find_angle(pos1, pos2=(0,0)):
    '''Expects one or two tuples or arrays with 2 as first shape, returns
    angle between them'''
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    a1 = np.arctan(pos1[1]/ pos1[0])
    a2 = np.arctan(pos2[1]/ pos2[0])
    return a1 - a2

def find_functions(system, planets = '', xy = False):
    if not planets:
        indexes = range(self.numberOfPlanets)
    else:
        indexes = planets.split(',')
        indexes = [int(i) for i in indexes]
    planetPos = np.load('positionsHomePlanet.npy')[:,indexes]
    times = np.load('times.npy')
    print planetPos.shape

    # angleFunction returns positions as function of theta, while
    # posFunction returns positions as function of time.
    posFunction = inter.interp1d(times[:-1], planetPos)
    theta, orbits_r = system.analytical_pos(planets)
    angleFunction = inter.interp1d(theta, orbits_r.T)
    return posFunction, angleFunction



def find_hohmann_params(system, A=0, B=1):
    '''returns parameters of the hohmann transfer from planet A to planet B
    at time t0 after start of sim (epoch)

    If only angle, it returns the required angle between the planets for a
    hohmann transfer. Else, currently it only returns the transfer time of
    the orbit
    
    :system: : Expects instance of class MySolarSystem
    : A, B : : Expects B > A
    '''
    pi = np.pi
    posFunction, angleFunction = find_functions(system, "%d, %d" %(A,B))
    n = 10
    r = np.zeros((n,2))
    for i,thet in enumerate(np.linspace(0,2*pi, n)):
        r[i] = angleFunction(thet)

    rA = np.mean(r, axis = 0)[0]
    rB = np.mean(r, axis = 0)[1]

    aTransfer = (rA+rB)/2.
    G = 4*pi**2
    M = system.starMass
    v_init_A        = np.sqrt(G*M/rA)
    v_final_B       = np.sqrt(G*M/rB)
    v_transfer_A    = np.sqrt(G*M*(2./rA - 1./aTransfer))
    v_transfer_B    = np.sqrt(G*M*(2./rB - 1./aTransfer))
    delta_VA        = v_transfer_A - v_init_A
    delta_VB        = v_final_B - v_transfer_B
    total_delta_V   = delta_VA + delta_VB
    e               = 1 -     rA/aTransfer #eccentrisity of transfer ellipse
    transferTime    = pi*np.sqrt((rA+rB)**3/(8*G*M))

    return {"v_init_A": v_init_A, 
            "v_final_B": v_final_B, 
            "total_delta_V":total_delta_V, 
            "delta_VA": delta_VA,
            "delta_VB": delta_VB,
            "e":e, 
            "transferTime":transferTime,
            "v_transfer_A":v_transfer_A,
            "v_transfer_B":v_transfer_B}
    
    
def angles_from_time(A,B,posFunction, angleFunction, t0):
    posA = posFunction(t0).T[A]
    rA = np.linalg.norm(posA)

    angleB = (find_angle(posA) + pi)%(2*pi)# The angle of pos of planet B
    posB = angleFunction(angleB % 2*pi)
    rB = np.linalg.norm(posB)
    
    transferTime = pi*np.sqrt((rA+rB)**3/(8*G*M))
    posB_at_final = posFunction(t0+transferTime)[B]
    angleB_at_final = find_angle(posB_at_final)
    return angleB, angleB_at_final

if __name__ == "__main__":
    import MySolarSystem as mySS
    system = mySS.MySolarSystem(87464)

    A = 0; B = 1
    times = np.linspace(0, system.period[A] + system.period[B], 100)
    posFunction, angleFunction = find_functions(system, A, B)
    aB, aBf= angles_from_time(A,B,posFunction,angleFunction,
                                        times)
    plt.plot(times,  aBf)
    plt.show()
    print aB, aBf

