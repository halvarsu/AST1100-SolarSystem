import numpy as np
import MySolarSystem as mss
system = mss.MySolarSystem(87464)


super_array = np.zeros((360,480,640,3))

phi = np.linspace(0,2*np.pi, 360)

for i,phi_0 in enumerate(phi):
    print 'Calculating deg: ',i, ', phi :', phi_0
    super_array[i] = system.projection(phi_0)

np.save('all_the_projections', super_array)
