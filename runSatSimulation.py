from MySateliteSim import MySateliteSim 

s = MySateliteSim(87464)
start_time = 14.25 #years
stop_time = start_time + 3. #years
init_boost = 9587#m/s
angle = 1.67*np.pi #radians
import time
sim_start = time.time()
dt_close = 1/(365.25*24*3600)
for x in range(3):
    s.launch_init(t0 = start_time, 
            init_boost = init_boost-x,init_angle = angle,
            ground_launch = True)
    pos, vel, acc, times= s.satelite_sim(
        t0=start_time, init_boost=init_boost-x, tN = stop_time, 
        init_angle = angle, dt_close = dt_close
    )
    print "length of pos ", pos.shape[0]
    print "Time used on sim: ", time.time() - sim_start,"s"
    norm = np.linalg.norm
    print "Last value: ", norm(pos[-1])
    print "Distance covered: ", norm(pos[0]-pos[-1])
    #s.plot()
    filenames = ['%s%d' %(name,init_vel-x) for name in ('vel','pos','times')]
    s.saveData(filenames)
