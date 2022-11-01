from gamma_distribution import Gamma_Variate_Function_Vector
import numpy as np
import matplotlib.pyplot as plt

class Gamma:
    def __init__(self, tmax, ymax, alpha, AT):
        self.tmax = tmax # x of the peak
        self.ymax = ymax # y of the peak
        self.alpha = alpha # shape parameter
        self.AT = AT # x shift

    def plot_gamma(self):
        t = np.arange(0, 1000, 0.1)
        plt.plot(t, Gamma_Variate_Function_Vector(t, self.tmax, self.ymax, self.alpha, self.AT), 'g')
        plt.show()

gamma1 = Gamma(50, 1, 1.0, 10)
gamma1.plot_gamma()