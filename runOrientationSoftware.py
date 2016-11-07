import MySateliteSim as m
s = m.MySateliteSim(87464)
from doppler import get_velocity

print s.orient('find_orient.png')
dl1 = float(raw_input('dl1: ') + 'e-9')
dl2 = float(raw_input('dl2: ') + 'e-9')
print get_velocity(dl1, dl2)
