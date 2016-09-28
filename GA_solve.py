from lambert_solve import init_planet, lambert_solve, lambert_plot
from MySateliteSim import MySateliteSim
import numpy as np

system = MySateliteSim(87464)
from PyKEP import AU, epoch

#planet1 = init_planet(0)
#planet4 = init_planet(4)

for x in range(7):
    planet = init_planet(x)

    eph = np.array(planet.eph(epoch(0)))
    print eph/AU
    print "-------------------------"
    print  x
    print "from pykep  x %7.4f  y %7.4f  " %(eph[0,0]/AU, eph[0,1]/AU)
    print "from module x %7.4f  y %7.4f  " %(system.x0[x], system.y0[x])


t1 =16.9105290335 
t2 =20.6651033379
A = 4
B = 5

l = lambert_solve(t1,t2,A,B)

#lambert_plot(t1,t2,A,B)
print l.get_v1()

