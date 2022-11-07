"""
TODO:  1) Check derived probabilities, 2) Reduction, 3) Final experimant with resulting plots
"""
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from gauss import Gauss
from mult_gauss import MultiGauss
from scipy.signal import find_peaks
import math

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

"""
ctr - mu
amp - p * 1/wid * 1/sqrt(pi)
wid - sqrt (2) * sigma
Note that if mean is not in x, the peak will not be properly visualised
"""
def gauss_func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        #print("i:")
        #print(i)
        #print(params)
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp(-((x - ctr)/wid)**2)
    return y

def build_multi_gauss_from_params(*params):
    m = MultiGauss([], [])
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        mean = ctr
        deviation = wid / math.sqrt(2) 
        g = Gauss(mean=mean, deviation=deviation)
        p = amp * wid * math.sqrt(math.pi)
        m.probabilities.append(p)
        m.gaussians.append(g)
    return m

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

def find_peaks_lib(x,y,min_height,width):
    peaks_x = []
    if width == None:
        peaks_i, _ = find_peaks(y, height=[min_height, 1])
    if min_height == None:
        peaks_i, _ = find_peaks(y, width=width)
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
    min_height = 0.0001
    peaks = find_peaks_lib(x,y,min_height=min_height,width=None)

    # when peaks are narrow --> calculate mean value
    if peaks == []:
        #peaks = [max_peak(x, y)]
        #peaks = find_peaks_lib(x,y,1)
        #peaks = filter_peaks(peaks, x, y, 50)
        mean = 0
        for i in range(len(x)):
            mean += x[i] * y[i]
        fit = gauss_func(x, *(mean, 1/math.sqrt(2*math.pi), math.sqrt(2)))
        m = build_multi_gauss_from_params([1.0], [Gauss(mean, 1)])
    else:
        print('peaks:')
        print(peaks)
        init_params = prepare_init_param(peaks)
        print(init_params)
        popt, pcov = curve_fit(gauss_func, x, y, p0 = init_params, maxfev=5000000)
        print('COV')
        print(pcov)
        print('parameters:')
        print(popt)
        fit = gauss_func(x, *popt)
        m = build_multi_gauss_from_params(*popt)
        # Avoiding flat fitting curves
        while max(fit) < 0.001:
            min_height *= 2
            peaks = find_peaks_lib(x,y,min_height=min_height,width=None)
            init_params = prepare_init_param(peaks)
            print(init_params)
            popt, pcov = curve_fit(gauss_func, x, y, p0 = init_params, maxfev=5000000)
            print('COV')
            print(pcov)
            print('parameters:')
            print(popt)
            fit = gauss_func(x, *popt)
            m = build_multi_gauss_from_params(*popt)
    #plt.plot(x, y)
    #plt.plot(x, fit , 'r-')
    m.plot_mult_gauss(x)
    print("Multi Gauss:")
    for i in range(len(m.gaussians)):
        print(m.gaussians[i].mean)
        print(m.probabilities[i])
    plt.show()

def max_value(x, y):
    max = 0
    for i in range(len(x)):
        if y(x) > max:
            max = y(x)
    return max
#mult_gauss1 = MultiGauss([0.5,0.5],[Gauss(10,2),Gauss(1,3)])
#t = np.arange(0, 20, 0.1)

#fit_gauss(t, mult_gauss1.mult_gauss_values())

def filter_peaks(peaks, x, y, max_number):
    filtered_peaks = []
    ind = collect_indexes_of_peaks(peaks, x)
    print("Indeces of peaks:")
    print(len(ind))
    for i in ind:
        cnt = 0
        for j in ind:
            if y[i] < y[j]:
                cnt += 1
        if cnt < max_number:
            filtered_peaks.append(x[i])
    print("Peaks after filter:")
    print(len(filtered_peaks))
    return filtered_peaks

def collect_indexes_of_peaks(peaks, x):
    ind = []
    for i in range(len(peaks)):
        for j in range(len(x)):
            if x[j]==peaks[i]:
                ind.append(j)
    return ind

def max_peak(x,y):
    max_i = 0
    max = 0
    for i in range(len(x)):
        if y[i] > max:
            max_i = i
            max = y[i]
    return x[max_i]

