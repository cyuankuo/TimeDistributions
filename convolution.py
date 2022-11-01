"""

@author: akalenkova
"""

#import numpy as np
#from math import exp
#import matplotlib.pyplot as plt 
#from scipy.optimize import minimize
#import ruptures as rpt
#from scipy import stats
from gauss import Gauss
from mult_gauss import MultiGauss

"""
Discere convolution of two functions f1 and f2 represented as lists
"""

def discrete_convolution(f1, f2):
    conv = []
    for i in range(len(f1)):
        for j in range(len(f2)):
            conv.append(f1[i]+f2[j])
    return conv

"""
Convolution of two gaussians
"""

def gauss_convolution(g1, g2):
    conv = Gauss(g1.mean+g2.mean, g1.variance+g2.variance)
    return conv

"""
Convolution of two sums of gaussians
mult1, mult2 are weighted sums of gaussians
returns another mult
"""

def mult_gauss_convolution(mult1, mult2):
    mult = MultiGauss([],[])
    for i in range(len(mult1.probabilities)):
        for j in range(len(mult2.probabilities)):
            mult.probabilities.append(mult1.probabilities[i]*mult2.probabilities[j])
            mult.gaussians.append(gauss_convolution(mult1.gaussians[i],mult2.gaussians[j]))
    return mult


gauss = mult_gauss_convolution(MultiGauss([0.0,1.0],[Gauss(10,2),Gauss(1,3)]), MultiGauss([0.0,0.1,0.9],[Gauss(1,2),Gauss(4,3),Gauss(7,15)]))
gauss.plot_mult_gauss()