from MySateliteSim import MySateliteSim
import numpy as np
import matplotlib.pyplot as plt
s = MySateliteSim(87464)

p = np.array((2,3,0))
v = np.array((1,-2,0))
R = 1
h = s.horizonVec(p,v,R)

t = np.linspace(0,2*np.pi,100)
line1 = np.array((p[:2],p[:2]+v[:2])).T
line2 = np.array((p[:2],p[:2]+h[:2])).T
line3 = np.array(((0,0),p[:2])).T
plt.plot(np.cos(t),np.sin(t))
plt.plot(line3[0],line3[1])
plt.plot(line1[0],line1[1])
plt.plot(line2[0],line2[1])
plt.axis('equal')
plt.show()

