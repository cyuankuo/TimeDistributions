import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

class Gauss:
    def __init__(self, mean, variance):
        self.mean = mean 
        self.variance = variance 

    def plot_gauss(self):
        t = np.arange(0, 20, 0.1)
        plt.plot(t, stats.norm.pdf(t, self.mean, self.variance))
        plt.show()

#gauss1 = Gauss(15,2)
#gauss1.plot_gauss()