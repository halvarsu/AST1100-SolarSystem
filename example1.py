import MySateliteSim as m
import numpy as np

system = m.MySateliteSim(87464)

launch = True
if launch:
    system.launch_init()
else:
    system.boost_init(pos0 = (3,3), vel0 = (-3,3), t0 = 0)

system.satelite_sim(speed_factor = [10,10,10])
system.plot()
