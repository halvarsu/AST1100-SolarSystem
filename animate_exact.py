# coding: utf-8
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig = plt.figure()
ax = fig.gca()

pos = np.load('positionsHomePlanet.npy')
positions = pos.swapaxes(1,2)
ax.plot(positions[0], positions[1], linewidth = 0.2)
ax.set_xlim([-28,28])
ax.set_ylim([-28,28])



length = pos.shape[-1]
res = 500

print length
print positions.shape

scat = plt.scatter(positions[0,0],positions[1,0], s = 50)

ax.axis('equal')

#ax.set_axis_bgcolor('k')#plt.set_background_color()
#sns.set_style('darkgrid', {'grid.color':0.5})
#ax.grid(False)

def update(i):
    j = (i*res)%length
    #print positions[:,j].shape, j

    scat.set_offsets(positions[:,j].T)

ani = FuncAnimation(fig, update, interval=5)
plt.show()
