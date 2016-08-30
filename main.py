"""
File: main.py
Author: Halvard Sutterud
Email: halvard.sutterud@gmail.com
Github: https://github.com/halvarsu
Description: The AST1100 project
"""

#from AST1100SolarSystem import  AST1100SolarSystem
from seed_reader import read_seed

    
class StarSystem(nil):
    """This is a doc string for this class"""
    def __init__(self, seed):
        super(StarSystem, self).__init__()
        self.seedl = seedl
        
def main():
    seed = read_seed('seed.dat')
    #myStarSystem = AST1100SolarSystem(seed)
    #starMass = myStarSystem.starMass
    #starRadius = myStarSystem.starRadius
    N = 3#myStarSystem.numberOfPlanets
    for i in xrange(N):
        #print 
    print seed
    return 0

if __name__ == "__main__":
    main()
