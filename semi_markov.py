import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special
from mult_gauss import MultiGauss
import math

class SemiMarkov:
    def __init__(self, states, transitions):
        self.states = states 
        # Each transition is a tuple: (from, to, probability, multi Gauss)
        self.transitions = transitions 

    