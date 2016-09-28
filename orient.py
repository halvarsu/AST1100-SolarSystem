"""
File: orient.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Given distance to planets with known positions, finds position
of satelite
"""

import numpy as np
time = 0
norm = np.linalg.norm
n = 2
#planetPos = np.array(((10.,6.),(-2.,5.)))#
planetPos = np.random.random((n,2)) * 10 - 5
pos = np.array((4.,5.))
#pos = np.random.random((2)) * 10 - 5
plt.scatter(pos[0],pos[1], c='y',s=100)
r_s = norm(pos)
theta = np.linspace(0,2*np.pi, 10000)
circle = r_s* np.array((np.sin(theta), np.cos(theta)))
plt.plot(circle[0], circle[1],color = 'k')


for i in range(n):
    x0 = planetPos[i,0]
    y0 = planetPos[i,1]
    r_p = norm(pos - planetPos[i])
    print x0, y0

    A = (x0**2 + y0**2 + r_s**2 - r_p**2)/(2*x0)
    b = y0/float(x0)
    C = np.sqrt(r_s**2*(b**2+1) - A**2)

    y = [(A*b + C)/(b**2+1),(A*b - C)/(b**2+1)]
    x = [A- y_i*b for y_i in y]
    for i in range(2):
        plt.scatter(x[i],
                y[i])
        print "comp", x[i], y[i] 
        circle = np.array((r_p*np.sin(theta)+x0, r_p*np.cos(theta)+y0))
        plt.plot(circle[0], circle[1]) 
        plt.plot((x0,pos[0]),(y0,pos[1]), '--')

        plt.axis('equal')
        plt.show()



