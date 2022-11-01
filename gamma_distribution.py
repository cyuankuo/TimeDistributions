# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 15:11:38 2019

@author: ancollet
"""

import numpy as np
from math import exp
import matplotlib.pyplot as plt 
from scipy.optimize import minimize
import ruptures as rpt
from scipy import stats

y_raw = []

"""
Build a vector of functio values
"""
def Gamma_Variate_Function_Vector(t, tmax, ymax, alpha, AT):
    f = []
    for i in range(len(t)):
        f.append(Gamma_Variate_Function_Madsen(t[i], tmax, ymax, alpha, AT))
    return f


def Gamma_Variate_Function_Madsen(t, tmax, ymax, alpha, AT):

    """
    Madsen, M. T., “A simplified formulation of the gamma variate function,
    ”Physics in Medicine and Biol-ogy37(7), 1597–1600 (1992).    

    "The fit of the gamma variate function has been used in numerous studies. 
    The main benefitsof the gamma variate function are the convenient 
    mathematical properties.  
    We found that the fitting of the gamma variate function leads to problems 
    because of the slow ”wash out” of the bolus".

    t: time value
    tmax: time of the peak
    ymax: peak intensity (absolute value)
    alpha: shape parameter
    AT: appearance time --> should be lower that tmax

    """
    if t <= AT:
        return 0
    if tmax <= AT:
        return float('nan')

    #t should always be positive
    t = abs((t - AT) / (tmax - AT))
    #print(t)
    f = ymax * pow(t, alpha) * exp(alpha * (1 - t))

    return f

def Multiple_Gamma_Variate_Function_Madsen(t, params):

    """
    Madsen, M. T., “A simplified formulation of the gamma variate function,
    ”Physics in Medicine and Biol-ogy37(7), 1597–1600 (1992).    

    t: time value
    tmax: time of the peak
    ymax: peak intensity (absolute value)
    alpha: shape parameter
    AT: appearance time --> should be lower that tmax

    """

    f = 0

    for i in range(int(len(params)/4)):
        tmax = params[i*4] 
        ymax = params[i*4 + 1]
        alpha = params[i*4 + 2]
        AT = params[i*4 + 3]

        if t <= AT:
            f += 0
        elif tmax <= AT:
            f += 0
        else:
            #t should always be positive
            t2 = abs((t - AT) / (tmax - AT))
            #print(t)
            f += ymax  * pow(t2, alpha) * exp(alpha * (1 - t2))
    return f
def Multiple_Gamma_Variate_Function_Madsen(t, params):

    """
    Madsen, M. T., “A simplified formulation of the gamma variate function,
    ”Physics in Medicine and Biol-ogy37(7), 1597–1600 (1992).    

    t: time value
    tmax: time of the peak
    ymax: peak intensity (absolute value)
    alpha: shape parameter
    AT: appearance time --> should be lower that tmax

    """

    f = 0

    for i in range(int(len(params)/4)):
        tmax = params[i*4] 
        ymax = params[i*4 + 1]
        alpha = params[i*4 + 2]
        AT = params[i*4 + 3]

        if t <= AT:
            f += 0
        elif tmax <= AT:
            f += 0
        else:
            #t should always be positive
            t2 = abs((t - AT) / (tmax - AT))
            #print(t)
            f += ymax  * pow(t2, alpha) * exp(alpha * (1 - t2))
    return f

# define objective function: SSE
def Objective_Fun(params, y_raw):

    # calculate y
    # calculate objective
    obj = 0.0
    for t in range(1, len(y_raw)):
        if y_raw[t] != 0: #and not isnan(y_raw[t]):
            y_fit = Multiple_Gamma_Variate_Function_Madsen(t, params)
            obj = obj + (y_fit-y_raw[t])**2
    # return result
    return obj

def Moving_Average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def Detect_Multiple_Gamma_Breakpoints(signal, n_bkps = 1, min_size = 10, 
                                      jump = 1, model = 'rbf',
                                      pen = 1):

    # Need to ensure a numpy array
    signal = np.array(signal)

    # Smoothing (to get rid of the noise), getting the base and removing it
    signal_smooth = Moving_Average(signal, n = min_size)
    print(signal_smooth)    

    # change point detection: Window segmentation
    algo = rpt.Window(width = min_size, model=model).fit(signal_smooth)

    # predict the ruptures with a given number of breakpoints
    breakpoints = algo.predict(n_bkps)

    # display the results
    rpt.display(signal, breakpoints)
    plt.plot(signal_smooth)
    plt.show()


    return breakpoints

def Find_Multiple_Gamma_Peaks_2(signal, n_peaks):

    # Case with only a single Gamma variate:
    if n_peaks == 1:
        i_max, j_max = 0, 0
        for i, j in enumerate(signal):
            if j > j_max:
                j_max = j
                i_max = i
        peaks_indices = {i_max: j_max}
        return peaks_indices

    # Else: Two Gamma to fit at least
    peaks_indices = {}

    #Detect the breakpoints in the series
    breakpoints = \
    Detect_Multiple_Gamma_Breakpoints(y_raw, n_bkps = n_peaks - 1, 
                                      min_size = 40, jump = 1, model = 'rbf',
                                      pen = 1)

    i1 = 0
    # Iterate over breakpoints
    for i2 in breakpoints:
        # Take the interval between breakpoints
        sub_signal = signal[i1 + 1 : i2 - 1]
        # Add the max of this sub_signal to the peak_index dico
        peak_index = \
        [i + i1 + 1 for i, j in enumerate(sub_signal) if j == max(sub_signal)]
        peaks_indices[peak_index[0]] = signal[peak_index[0]]
        # Set the lower bound for the next iteration as the current 
        # interval upperbound
        i1 = i2

    return peaks_indices

def Fit_TS_Multiple_Gamma_Variates(y_raw, n_peaks):

    # Find the peaks
    peaks = Find_Multiple_Gamma_Peaks_2(y_raw, n_peaks)
    # Initiate a preset of parameters from the peaks found
    params_raw = Initiate_Parameters(peaks)
    # cons = f(nb of peaks)
    cons = Initiate_Parameters_Constraints(peaks)
    # bonds = f(nb of peaks)
    bnds = Initiate_Parameters_Bounds(peaks)

    tol = 1e-7

    """
    #Solving
    solution = minimize(Objective_Fun, x0 = params, method='CG', bounds=bnds,
                        args = (y_raw), constraints = None)
    """

    """
    #Solving
    solution = minimize(Objective_Fun, x0 = params, method='SLSQP', bounds=bnds,
                        args = (y_raw), constraints = cons)

    """

    # Solving
    solution = minimize(Objective_Fun, x0 = params_raw, method='COBYLA', 
                        bounds=bnds, args = (y_raw), tol = tol)

    # Get the optimized parameters
    params_fit = solution.x

    # Generate the fitted curve (general + individual)
    y_fits = Get_Y_Fits(y_raw, params_fit, n_peaks)

    return y_fits, params_raw, params_fit

def Get_Y_Fits(y_raw, params, n_peaks):

    # Store the curves as a dico
    y_fits = {}

    y_fit = np.zeros(len(y_raw))
    for t in range(len(y_raw)):
        y_fit[t] = Multiple_Gamma_Variate_Function_Madsen(t, params)

    y_fits['Global gamma fit curve'] = y_fit

    if n_peaks > 1:
        for pk in range(1, n_peaks + 1):
            y_fit = np.zeros(len(y_raw))
            sub_params = params[(pk - 1) * 4: pk * 4] 
            #print(sub_params)
            for t in range(len(y_raw)):
                y_fit[t] = \
                Multiple_Gamma_Variate_Function_Madsen(t, sub_params)
            y_fits['Gamma fit cruve ' + str(pk)] = y_fit

    return y_fits


def Initiate_Parameters(peaks):

    # Initiate the parameters based on the found curves
    params = np.zeros(4 * len(peaks))    
    for pk in range(len(peaks)):         
        # initial guesses
        params[pk * 4 + 0] = list(peaks)[pk]                    # tmax
        params[pk * 4 + 1] = peaks[params[pk * 4 + 0]] * 0.7    # ymax
        params[pk * 4 + 2] = 1                                  # alpha
        params[pk * 4 + 3] = params[pk * 4 + 0] - 1000          # AT

    return params

def Initiate_Parameters_Constraints(peaks):

    #Empty list
    cons = []
    #Setting the constraints as dictionnaries
    for pk in range(len(peaks)):
        cons.append({'type': 'ineq', 
                     'fun': lambda x: x[pk * 4 + 0] - x[pk * 4 + 3]})
    return cons

def Initiate_Parameters_Bounds(peaks):

    bnds = ()
    for pk in range(len(peaks)):        
        # bounds on variables
        # tmax belongs to R
        # ymax belongs to R+
        # alpha belongs to R+
        # AT belongs to R
        bnds += ((-10000, 10000), (0, 5000), (0, 50), (-10000, 10000))

    return bnds

def Get_Variance_Explaination_R2(y_raw, y_fit):

    slope, intercept, r_value, p_value, std_err = \
    stats.linregress(y_raw, y_fit)
    r2 = r_value**2

    return r2

def Display(r2, params_raw, params_fit):

    # show initial objective
    #print('  • Initial Objective score: ' \
    #      + str(Objective_Fun(params_raw, y_raw)))

    # show final objective
    #print('  • Final Objective score:   ' \
    #      + str(Objective_Fun(params_fit, y_raw)))

    cR2 = "  • R² correlation = " + str(r2)

    print(cR2)

    n_peaks = int(len(params_fit)/4)

    for pk in range(n_peaks):
        # print solution
        print('Solution for gamma cuvrve', str(pk + 1), ':')
        ctmax =  '  • tmax   = ' + str(params_fit[pk * 4])
        print(ctmax)
        cymax =  '  • ymax   = ' + str(params_fit[pk * 4 + 1])
        print(cymax)
        calpha = '  • alpha  = ' + str(params_fit[pk * 4 + 2])
        print(calpha)
        cAT =    '  • AT     = ' + str(params_fit[pk * 4 + 3])
        print(cAT)

def fit_gamma(y):
    y_raw = []
    x_raw = range(1000)
    for i in x_raw:
        if i in y.keys():
            y_raw.append(y.get(i))
        else:
            y_raw.append(0)
    #------
    print("Y RAW:")
    print(y_raw)
    plt.plot(y_raw)
    plt.show() 
    #------
    for n_peaks in [1, 2, 3, 4, 5]:

        print('Fitting the curve with ', n_peaks, ' Gamma functions...')

        #Fit the functions
        y_fits, params_raw, params_fit = \
        Fit_TS_Multiple_Gamma_Variates(y_raw, n_peaks)

        #Calculate R2
        r2 = \
        Get_Variance_Explaination_R2(y_raw, y_fits['Global gamma fit curve'])

        #Display the regression parameters and info:
        Display(r2, params_raw, params_fit)

        #Plot the results
        plt.plot(y_raw)        
        for key in y_fits.keys():
            plt.plot(y_fits[key])
        plt.show()

#Test
if __name__ == '__main__':

    # Generate a multi - Gamma variate curve
    x_raw = range(1000)

    #params = [100, 450, 0.9, -100, 250, 300, 1.2, 100, 600, 300, 2, 400]

    # y_raw = [Multiple_Gamma_Variate_Function_Madsen(x, params) + \
    #            float(np.random.normal(0,2,1)) for x in x_raw]
    #y_raw = [Multiple_Gamma_Variate_Function_Madsen(x, params)  for x in x_raw]
    #y_raw = [Multiple_Gamma_Variate_Function_Madsen(x, params) for x in x_raw]
    y_raw = [float(np.random.normal(10,2,1)) for x in x_raw]
    #Show the raw curve
    print('raw signal:')
    plt.plot(y_raw) #.fig_title = 'Raw noisy data'
    plt.show()

    # Fit successively with one, two and then three Gamma - Variate
    for n_peaks in [1, 2, 3, 4, 5]:

        print('Fitting the curve with ', n_peaks, ' Gamma functions...')

        #Fit the functions
        y_fits, params_raw, params_fit = \
        Fit_TS_Multiple_Gamma_Variates(y_raw, n_peaks)

        #Calculate R2
        r2 = \
        Get_Variance_Explaination_R2(y_raw, y_fits['Global gamma fit curve'])

        #Display the regression parameters and info:
        Display(r2, params_raw, params_fit)

        #Plot the results
        plt.plot(y_raw)        
        for key in y_fits.keys():
            plt.plot(y_fits[key])
        plt.show()
