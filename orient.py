"""
File: orient.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: Given distance to planets with known positions, finds position
of satelite
"""

import matplotlib.pyplot as plt
import numpy as np
import MySateliteSim as m
if __name__ == "__main__":
    system = m.MySateliteSim(87464)
    norm = np.linalg.norm
    posFunc = system.getPositionsFunction()


    #time = 4.784784784
    #time = 0.964020964021
    #time = 10.499049905
    planetPos = posFunc(2).T
    t = np.linspace(0,20, 100000)
    pP = posFunc(t)[:,0]
    print pP.shape
    pN = norm(pP, axis = 0)
    p_list = np.load('pos.npy')
    print np.amin(pN- p_list[-1])
    print np.amax(pN)

    
    arg =  np.argmin(np.abs(pN- p_list[-1]))
    sort =np.argsort(np.abs(pN- p_list[-1])) 
    print t[arg]
    best_times =  t[sort][:10]
    print best_times


    #import sys
    #sys.exit()
    pos = np.random.random((2)) * 20 - 10
    relpos = planetPos - pos 
    rel_r = norm(relpos, axis = 1)
    #p_list = np.append(rel_r, norm(pos))
    p_list = np.load('pos.npy')

    time = 14.3199999998
    m.find_distance(p_list, time = time)

    print p_list
