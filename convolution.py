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
import math

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
    conv = Gauss(g1.mean+g2.mean, math.sqrt((g1.deviation)**2+(g2.deviation)**2))
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
    mult.remove_small_prob_gauss(0.001)
    return mult

def mult_gauss_self_convolution(mult1, k):
    mult = MultiGauss([],[])
    for i in range(k):
        mult = mult_gauss_convolution(mult, mult1)
    mult.remove_small_prob_gauss(0.001)
    return mult

def mult_gauss_sum(mult1, mult2, p1, p2):
    sum = MultiGauss([],[])
    for i in range(len(mult1.probabilities)):
        sum.probabilities.append(p1*mult1.probabilities[i])
        sum.gaussians.append(mult1.gaussians[i])
    
    for i in range(len(mult2.probabilities)):
        sum.probabilities.append(p2*mult2.probabilities[i])
        sum.gaussians.append(mult2.gaussians[i])
    
    sum.remove_small_prob_gauss(0.001)
    return sum


#gauss = mult_gauss_convolution(MultiGauss([0.0,1.0],[Gauss(10,2),Gauss(1,3)]), MultiGauss([0.0,0.1,0.9],[Gauss(1,2),Gauss(4,3),Gauss(7,15)]))
#gauss.plot_mult_gauss()