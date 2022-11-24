import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special
from gauss import Gauss
import math
from copy import deepcopy


class MultiGauss:
    def __init__(self, probabilities, gaussians):
        self.probabilities = probabilities 
        self.gaussians = gaussians 

    def plot_mult_gauss(self, x):
        f = [0] * len(x)
        for i in range(len(self.probabilities)):
            print(self.probabilities[i])
            print(self.gaussians[i].mean)
            print(self.gaussians[i].deviation)
            print(self.probabilities[i])
            f += self.probabilities[i] * stats.norm.pdf(x, self.gaussians[i].mean, self.gaussians[i].deviation)
            #plt.plot(x, f)
            #plt.show()
            #f = stats.norm.pdf(x, self.gaussians[i].mean, self.gaussians[i].deviation)
        plt.plot(x, f)
        #plt.show()
    def remove_out_bounds_gauss(self, x):
        i =0 
        while i < len(self.probabilities):
            if (self.gaussians[i].mean > max(x)):
                del self.probabilities[i]
                del self.gaussians[i]
            else:
                i+=1
        self.normalise_gauss()

    def mult_gauss_values(self):
        t = np.arange(0, 20, 0.1)
        f = [0] * len(t)
        for i in range(len(self.probabilities)):
            f += self.probabilities[i] * stats.norm.pdf(t, self.gaussians[i].mean, self.gaussians[i].deviation)
        return f
    
    def normalise_gauss(self):
        sum_p = 0
        for i in range(len(self.probabilities)):
            sum_p += self.probabilities[i]

        for i in range(len(self.probabilities)):
            self.probabilities[i] /= sum_p

        return self

    def truncate_gauss(self, threshold):
        for i in range(len(self.gaussians)):
            negative_area = self.gaussians[i].calc_negative_area()
            print(negative_area)
            if  negative_area > threshold:
               self.gaussians[i].deviation = - (self.gaussians[i].mean/(math.sqrt(2)*special.erfinv(2*threshold - 1)))
        return self
    
    def remove_small_prob_gauss(self, threshold):
        i =0
        while i < len(self.probabilities):
            if (self.probabilities[i] < threshold):
                del self.probabilities[i]
                del self.gaussians[i]
            else: 
                i += 1
        self.normalise_gauss()
    
    def remove_zero(self):
        i=0
        while i < len(self.probabilities):
            if (self.gaussians[i].mean == 0) and (self.gaussians[i].deviation == 0):
                del self.probabilities[i]
                del self.gaussians[i]
            else:
                i += 1
        self.normalise_gauss()


#mult_gauss1 = MultiGauss([0.3,0.5],[Gauss(10,2),Gauss(1,3)])
#mult_gauss1.plot_mult_gauss()
#mult_gauss1 = mult_gauss1.normalise_gauss()
#mult_gauss1.plot_mult_gauss()
#mult_gauss1.truncate_gauss(0.05)
#mult_gauss1.plot_mult_gauss()
#plt.show()