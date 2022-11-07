import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special
from gauss import Gauss
import math

class MultiGauss:
    def __init__(self, probabilities, gaussians):
        self.probabilities = probabilities 
        self.gaussians = gaussians 

    def plot_mult_gauss(self):
        t = np.arange(0, 20, 0.1)
        f = [0] * len(t)
        for i in range(len(self.probabilities)):
            f += self.probabilities[i] * stats.norm.pdf(t, self.gaussians[i].mean, self.gaussians[i].deviation)
            #print(f)
        plt.plot(t, f)

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

mult_gauss1 = MultiGauss([0.3,0.5],[Gauss(10,2),Gauss(1,3)])
mult_gauss1.plot_mult_gauss()
mult_gauss1 = mult_gauss1.normalise_gauss()
mult_gauss1.plot_mult_gauss()
mult_gauss1.truncate_gauss(0.05)
mult_gauss1.plot_mult_gauss()
plt.show()