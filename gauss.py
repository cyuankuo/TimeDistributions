import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math

class Gauss:
    def __init__(self, mean, deviation):
        self.mean = mean 
        self.deviation = deviation 

    def plot_gauss(self):
        t = np.arange(0, 20, 0.1)
        #plt.plot(t, stats.norm.pdf(t, self.mean, self.deviation))
        #plt.show()

    """
    Calculate the areaunder the curve when x < 0
    """
    def calc_negative_area(self):
        cdf = 0.5 * (1 + math.erf(-self.mean/(self.deviation*math.sqrt(2))))
        return cdf
#gauss1 = Gauss(15,2)
#gauss1.plot_gauss()