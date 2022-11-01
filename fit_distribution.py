"""
TODO: 1) Finish with the built in method, 2) Truncate Gauss, 3) Reduction, 4) Final experimants for 10 logs
"""
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from gauss import Gauss
from mult_gauss import MultiGauss
from scipy.signal import find_peaks

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

def find_peaks_custom(x,y,max_number):
    peaks = []
    max_x = len(x)
    print(max_x)
    print(max_x//max_number)
    for i in range (0, max_x - max_x//max_number, max_x//max_number):
        max = 0
        peak_x = 0
        for j in range(i,i + max_x//max_number):
            if y[j] > max:
                max = y[j]
                peak_x = x[j]
        peaks.append(peak_x)
    return peaks

def find_peaks_lib(x,y):
    peaks_x = []
    peaks_i, _ = find_peaks(y,width=20)
    print('Peaks I:')
    print(len(peaks_i))
    for peak_i in peaks_i:
        for i in range(len(x)):
            if (i==peak_i):
                peaks_x.append(x[i])
    return peaks_x

def prepare_init_param(peaks):
    init_params = []
    for peak in peaks:
        init_params.append(peak)
        init_params.append(1)
        init_params.append(1)
    return init_params

def fit_gauss(x,y):
    print('x:')
    print(x)
    print('y:')
    print(y)
    #y = moving_average(y, 3)
    #print(y)
    peaks = find_peaks_lib(x,y)
    print('peaks:')
    print(peaks)
    init_params = prepare_init_param(peaks)
    print(init_params)
    popt, pcov = curve_fit(gauss_func, x, y, p0 = init_params, maxfev=5000000)
    print('parameters:')
    print(popt)
    fit = gauss_func(x, *popt)

    plt.plot(x, y)
    plt.plot(x, fit , 'r-')
    plt.show()

#mult_gauss1 = MultiGauss([0.5,0.5],[Gauss(10,2),Gauss(1,3)])
#t = np.arange(0, 20, 0.1)

#fit_gauss(t, mult_gauss1.mult_gauss_values())

