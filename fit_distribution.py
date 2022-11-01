from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from gauss import Gauss
from mult_gauss import MultiGauss

"""
Smooth a function
"""
def moving_average(y, n=3) :
    sum = np.zeros(len(y))
    for i in range(len(y)-n+1):
        sum[i] = 0
        for j in range (i,i+n-1):
            sum[i] += y[j]
        sum[i] /= 3
    y = sum
    return y 

def gauss_func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

def fit_gauss(x,y):
    print(x)
    print(y)
    y = moving_average(y, 3)
    print(y)
    popt, pcov = curve_fit(gauss_func, x, y, p0 = [1,1,1,10,1,1], maxfev=1000000)
    print(popt)
    fit = gauss_func(x, *popt)

    plt.plot(x, y)
    plt.plot(x, fit , 'r-')
    plt.show()

#mult_gauss1 = MultiGauss([0.5,0.5],[Gauss(10,2),Gauss(1,3)])
#t = np.arange(0, 20, 0.1)

#fit_gauss(t, mult_gauss1.mult_gauss_values())

