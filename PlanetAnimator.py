"""
File: animate_exact.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Arguments should be indexes of wanted planets
"""
    
# coding: utf-8
#import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from MySolarSystem import MySolarSystem
from matplotlib.animation import FuncAnimation
import sys
from time import sleep

try:
    if len(sys.argv) ==1:
        raise IndexError
    elif len(sys.argv) > 8:
        raise ValueError
    planets = []
    for p in sys.argv[1:]:
        planets.append(int(p))
except IndexError:
    planets = range(7)
except ValueError:
    print ("bad usage, expects int arguments (max 7)(planet indexes)")
    sys.exit()


class PlanetAnimator(MySolarSystem):

    """Docstring for PlanetAnimator. """

    def __init__(self, seed, years=0):
        MySolarSystem.__init__(self, seed)
        """TODO: to be defined1."""
        self.years = years
        self.sleep_time = 0
        self.forwards = True
        self.j = 0
        self.res = 200
        self.direction= "Forwards"
        self.paused = False
        

    def load_files(self, times = 'times.npy', sat=False):
        pos = np.load('positionsHomePlanet.npy')
        self.times = np.load('times.npy')
        self.positions = pos[:,planets].swapaxes(0,1).swapaxes(1,2)
        self.plot_pos = np.load('analyticalPositions.npy')
        self.length = pos.shape[-1]


    def animate(self, save = False):
        self.fig = plt.figure(figsize = (7,7))
        self.ax = self.fig.gca()
        self.ax.plot(self.plot_pos[0].T, self.plot_pos[1].T, linewidth = 0.2)
        positions = self.positions
        color = np.random.random((3))
        self.scat = plt.scatter(positions[:,0,0],positions[:,0,1], s = 50,c=color)

        self.fig.suptitle("Sleep time: %f"%self.sleep_time)
        self.max_dist = np.max(np.abs(positions)) + 1
        self.ax.axis('equal')
        self.ax.set_xlim([-self.max_dist,self.max_dist])
        self.ax.set_ylim([-self.max_dist,self.max_dist])
        self.ani = FuncAnimation(self.fig, self.update, interval=5)
        cid = self.fig.canvas.mpl_connect('key_press_event', self.onpress)
        plt.show()
        if save:
            from JSAnimation import HTMLWriter
            self.ani.save('animation.html', writer=HTMLWriter(embed_frames=True))


    def update(self,i):
        if self.forwards:
            self.j += self.res
            if self.j >= self.length:
                self.j -= self.length
        else:
            self.j -= self.res
            if self.j < 0:
                self.j += self.length
        self.scat.set_offsets(self.positions[:,self.j])
        self.fig.suptitle("Time = %f, Direction = %s, speed = %d"%\
                ( self.times[self.j], self.direction, self.res))


    def onpress(self, event):
        key = event.key
        if key == " ":
            if self.paused:
                self.ani.event_source.start()
                self.paused = False
            else:
                self.ani.event_source.stop()
                self.paused = True
        elif key == "w":
            self.max_dist *= 0.5
            self.ax.set_xlim([-self.max_dist,self.max_dist])
            self.ax.set_ylim([-self.max_dist,self.max_dist])
        elif key == "e":
            self.max_dist *= 2
            self.ax.set_xlim([-self.max_dist,self.max_dist])
            self.ax.set_ylim([-self.max_dist,self.max_dist])
        elif key == "right":
            self.res += 50
        elif key == "left":
            self.res -= 50
        elif key == "down" or key == "up":
            self.forwards = not self.forwards
        elif key == "r":
            self.j = 0
        self.direction ="Forwards" if self.forwards else "Backwards"




if __name__ == "__main__":
    pa = PlanetAnimator(87464)
    pa.load_files()
    pa.animate()

    
