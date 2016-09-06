# coding: utf-8
from AST1100SolarSystem import AST1100SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class MySolarSystem(AST1100SolarSystem):

    """This is a doc string"""

    def __init__(self, seed):
        """TODO: to be defined1.

        :seed: TODO

        """
        AST1100SolarSystem.__init__(self, seed)

        self._seed = seed

    def polar_euler(self):
        for i in range(self.numberOfPlanets):
            print self.mass[i]

        
if __name__ == "__main__":
    seed = 87464
    system = MySolarSystem(seed)
    system.polar_euler()
