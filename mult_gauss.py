import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from gauss import Gauss

class MultiGauss:
    def __init__(self, probabilities, gaussians):
        self.probabilities = probabilities 
        self.gaussians = gaussians 

    def plot_mult_gauss(self):
        t = np.arange(0, 20, 0.1)
        f = [0] * len(t)
        for i in range(len(self.probabilities)):
            f += self.probabilities[i] * stats.norm.pdf(t,  self.gaussians[i].mean, self.gaussians[i].variance)
            print(f)
        plt.plot(t, f)
        plt.show()

    def mult_gauss_values(self):
        t = np.arange(0, 20, 0.1)
        f = [0] * len(t)
        for i in range(len(self.probabilities)):
            f += self.probabilities[i] * stats.norm.pdf(t,  self.gaussians[i].mean, self.gaussians[i].variance)
        return f

#mult_gauss1 = MultiGauss([0.5,0.5],[Gauss(10,2),Gauss(1,3)])
#mult_gauss1.plot_mult_gauss()