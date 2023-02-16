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

    def plot_mult_gauss(self, x, label, color):
        f = [0] * len(x)
        for i in range(len(self.probabilities)):
            #print(self.probabilities[i])
            f += self.probabilities[i] * stats.norm.pdf(x, self.gaussians[i].mean, self.gaussians[i].deviation)
            #plt.plot(x, f)
            #plt.show()
            #f = stats.norm.pdf(x, self.gaussians[i].mean, self.gaussians[i].deviation)
#        if label == 'print!':
#            plt.clf()
#            frame1 = plt.gca()
#            frame1.axes.xaxis.set_ticklabels([])
#            frame1.axes.yaxis.set_ticklabels([])
#            frame1.axes.xaxis.set_ticks([])
#            frame1.axes.yaxis.set_ticks([])
        plt.plot(x, f, color, linewidth=2, label=label)

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
            if  negative_area > threshold:
                #print(negative_area)
                self.gaussians[i].deviation = - (self.gaussians[i].mean/(math.sqrt(2)*special.erfinv(2*threshold - 1)))
        return self
    
    def remove_small_prob_gauss(self, threshold):
        #print("Old size:")
        #print(len(self.probabilities))
        i =0
        while i < len(self.probabilities):
            if (self.probabilities[i] < threshold):
                del self.probabilities[i]
                del self.gaussians[i]
            else: 
                i += 1
        self.normalise_gauss()
        #print("New size:")
        #print(len(self.probabilities)) 

    
    def unify_small_prob_gauss(self, threshold):
        #print("Old size:")
        #print(len(self.probabilities))
        i = 0
        sum_prob = 0
        sum_means = 0
        sum_deviations = 0
        while i < len(self.probabilities):
            if (self.probabilities[i] < threshold):
                sum_prob += self.probabilities[i]
            i += 1
        #print("Sum probabilitties:")
        #print(sum_prob)     
        if sum_prob > 0:
            i = 0
            while i < len(self.probabilities):
                if (self.probabilities[i] < threshold):
                    sum_means += self.probabilities[i]/sum_prob*self.gaussians[i].mean
                    sum_deviations += self.probabilities[i]/sum_prob*(self.gaussians[i].mean**2 + self.gaussians[i].deviation**2)    
                    del self.probabilities[i]
                    del self.gaussians[i]
                else:
                    i += 1
            
        sum_deviations = np.sqrt(sum_deviations - sum_means**2)
        if sum_prob > 0:
            self.probabilities.append(sum_prob)  
            self.gaussians.append(Gauss(sum_means, sum_deviations))
        #print("Sum probabilitties:")
        #print(sum_prob) 
        #print("New size:")
        #print(len(self.probabilities)) 
    
    def remove_zero(self):
        i=0
        while i < len(self.probabilities):
            if (self.gaussians[i].mean == 0) and (self.gaussians[i].deviation == 0):
                del self.probabilities[i]
                del self.gaussians[i]
            else:
                i += 1
        self.normalise_gauss()
    
    def calculate_mean(self):
        i =0
        mean = 0 
        while i < len(self.probabilities):
            mean += self.probabilities[i]*self.gaussians[i].mean
            i += 1
        return mean
    
    def calculate_mode(self):
        i =0
        mode = 0
        max_prob = 0 
        while i < len(self.probabilities):
            if max_prob < self.probabilities[i]:
                mode = self.gaussians[i].mean
                max_prob = self.probabilities[i]
            i += 1
        return mode
    
    def calculate_peak(self):
        i =0
        peak = 0
        deviation = 0
        prob = 0 
        while i < len(self.probabilities):
            prob = self.probabilities[i]
            peak = self.gaussians[i].mean
            deviation = self.gaussians[i].deviation
            if peak > 10 and peak < 20:
                break
            i += 1
        return peak, deviation, prob
    
    def calculate_sum_probabilities(self):
        i =0
        prob = 0 
        while i < len(self.probabilities):
            prob += self.probabilities[i]
            i += 1
        return prob

#mult_gauss1 = MultiGauss([0.3,0.5],[Gauss(10,2),Gauss(1,3)])
#mult_gauss1.plot_mult_gauss()
#mult_gauss1 = mult_gauss1.normalise_gauss()
#mult_gauss1.plot_mult_gauss()
#mult_gauss1.truncate_gauss(0.05)
#mult_gauss1.plot_mult_gauss()
#plt.show()