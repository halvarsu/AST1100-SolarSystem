import sys
import subprocess
from astropy import units as u

def main():
    filename = 'seed.dat'
    print read_seed(filename)
    username = 'halvarsu'
    print get_seed(username)
    


def read_seed(filename):
    infile = open(filename)
    line = infile.readline()
    return int(line)

def get_seed(username):
    output = subprocess.Popen(["python", "myseed.pyc", username], stdout=subprocess.PIPE).communicate()[0]
    return int(output)


if __name__ == "__main__":
    main()
