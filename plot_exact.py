# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

pos = np.load('positionsHomePlanet.npy')
positions = pos.swapaxes(1,2)
plt.plot(positions[0], positions[1])
plt.axis('equal')
plt.show()
