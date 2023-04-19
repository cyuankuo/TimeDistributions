import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special
from gauss import Gauss
import math
from copy import deepcopy
from scipy.stats import truncnorm


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
    
    def calculate_peaks(self):
        i =0
        prob = 0
        max_probabilities = []
        max_means = []
        max_length = 5 
        while i < len(self.probabilities):
            prob = self.probabilities[i]
            if i < max_length:
                max_probabilities.append(prob)
            else:
                for max_prob in max_probabilities:
                    if prob > max_prob:
                        max_probabilities.remove(max_prob)
                        max_probabilities.append(prob)
                        break
            i += 1
        i = 0
        while i < len(max_probabilities):
            max_means.append(self.gaussians[i].mean)
            i += 1

        return max_means, max_probabilities
    
    def calculate_sum_probabilities(self):
        i =0
        prob = 0 
        while i < len(self.probabilities):
            prob += self.probabilities[i]
            i += 1
        return prob
    
    def calc_observed(self, lower_bound, upper_bound, samp):
        cnt = 0
        for s in samp:
            if s <= upper_bound and s >= lower_bound:
                cnt += 1
        return cnt

    def calc_expected(self, lower_bound, upper_bound, samp):
        expected = 0
        for i in range(len(self.probabilities)):
            cdf_lower, cdf_upper = stats.norm.cdf([lower_bound, upper_bound], self.gaussians[i].mean, self.gaussians[i].deviation)
            expected += len(samp)*self.probabilities[i] * (cdf_upper-cdf_lower)
        return expected


    def calc_expected_truncated(self, lower_bound, upper_bound, samp):
        expected = 0
        for i in range(len(self.probabilities)):
            a, b = (0 - self.gaussians[i].mean) / self.gaussians[i].deviation, np.inf
            cdf_lower, cdf_upper = truncnorm.cdf([lower_bound, upper_bound], a=a, b=b, loc=self.gaussians[i].mean, scale=self.gaussians[i].deviation)
            expected += len(samp)*self.probabilities[i] * (cdf_upper-cdf_lower)
        return expected

    def calc_chi_square(self, bins, samp):
        chi_square = 0
        #print(samp)
        upper_bound = np.max(samp)
        step = np.max(samp)/bins

        for i in range(bins):
           lower_bound = i*step
           upper_bound = (i+1)*step
           #print(lower_bound)
           #print(upper_bound)
           observed = self.calc_observed(lower_bound, upper_bound, samp)
           expected = self.calc_expected_truncated(lower_bound, upper_bound, samp)
           #print()
           #print(observed)
           #print(expected)
           #print((observed-expected)**2 / expected)
           chi_square += (observed-expected)**2 / expected
           

        return chi_square

    def calc_kl_divergence(self, bins, samp):
        kl_divergence = 0
        #print(samp)
        upper_bound = np.max(samp)
        step = np.max(samp)/bins

        for i in range(bins):
           lower_bound = i*step
           upper_bound = (i+1)*step
           #print(lower_bound)
           #print(upper_bound)
           #P()
           observed = self.calc_observed(lower_bound, upper_bound, samp)/len(samp)
           expected = self.calc_expected_truncated(lower_bound, upper_bound, samp)/len(samp)
           #print()
           #print(observed)
           #print(expected)
           #print((observed-expected)**2 / expected)
           #print(expected)
           #print(observed)
           if observed != 0 and expected != 0:
                kl_divergence += observed*(np.log(observed/expected))
           

        return kl_divergence
    
    def calc_chi_square_uniform(self, bins, samp):
        chi_square = 0
        upper_bound = np.max(samp)
        step = np.max(samp)/bins

        for i in range(bins):
           lower_bound = i*step
           upper_bound = (i+1)*step
           #print(lower_bound)
           #print(upper_bound)
           observed = self.calc_observed(lower_bound, upper_bound, samp)
           expected = step
           #print()
           #print(observed)
           #print(expected)
           #print((observed-expected)**2 / expected)
           chi_square += (observed-expected)**2 / expected
           

        return chi_square


    def plot_trunc_mult_gauss(self, x, label, color):
        f = [0] * len(x)
        for i in range(len(self.probabilities)):
            #print(self.probabilities[i])
            #a, b = 0, np.inf
            a, b = (0 - self.gaussians[i].mean) / self.gaussians[i].deviation, np.inf
            f += self.probabilities[i] * truncnorm.pdf(x, a=a, b=b, loc=self.gaussians[i].mean, scale=self.gaussians[i].deviation)
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
#mult_gauss1 = MultiGauss([0.3,0.5],[Gauss(10,2),Gauss(1,3)])
#mult_gauss1.plot_mult_gauss()
#mult_gauss1 = mult_gauss1.normalise_gauss()
#mult_gauss1.plot_mult_gauss()
#mult_gauss1.truncate_gauss(0.05)
#mult_gauss1.plot_mult_gauss()
#plt.show()