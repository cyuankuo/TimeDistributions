from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.objects.dfg.exporter import exporter as dfg_exporter
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery
from pm4py.statistics.start_activities.log import get as start_activities
from pm4py.statistics.end_activities.log import get as end_activities
from pm4py.visualization.dfg import visualizer as dfg_visualization
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import statsmodels.api as sm
from itertools import groupby
import gamma_distribution
import numpy as np
import seaborn as sns
from fit_distribution import fit_gauss
import collections

def extract_times():
    for trace in log:
        first = True
        for next_event in trace:
            if not first:
                time = next_event['time:timestamp'] - event['time:timestamp']
                hours = time.total_seconds()//3600
                if hours > 250:
                    continue
                if not event['concept:name'] in times_dictionary.keys():
                    times_dictionary[event['concept:name']] = [time.total_seconds()//3600]
                else:
                    times_dictionary[event['concept:name']].append(time.total_seconds()//3600)
            event = next_event
            first = False

def extract_times_with_future():
    for trace in log:
        first = True
        for next_event in trace:
            if not first:
                time = next_event['time:timestamp'] - event['time:timestamp']
                if not event['concept:name'] + '+' + next_event['concept:name'] in times_dictionary.keys():
                    times_dictionary[event['concept:name'] + '+' + next_event['concept:name']] = [time.total_seconds()//3600]
                else:
                    times_dictionary[event['concept:name'] + '+' + next_event['concept:name']].append(time.total_seconds()//3600)
            event = next_event
            first = False
# retrieve distribution of values from given y
def retrieve_distribution(y):
    result = {}
    for i in y:
        count_i = y.count(i)
        result[i] = count_i
        result[i] /= len(y)
    #gamma_distribution.fit_gamma(result)
    return result

variant = xes_importer.Variants.ITERPARSE
parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
log = xes_importer.apply('/Users/a1230101//Documents/GitHub/TimeDistributions/logs/DomesticDeclarations.xes', 
    variant=variant, parameters=parameters)
times_dictionary = {}

extract_times_with_future()

"""
Fitting using Gaussian KDE
"""

fig = plt.figure(figsize=(10,220))
i = 1
for key in sorted(times_dictionary.keys()):
    #dist = retrieve_distribution(times_dictionary.get(key))
    h = fig.add_subplot(len(times_dictionary), 1, i)
    y, x, _ = h.hist(times_dictionary.get(key), bins=1000, facecolor='g', density = True,)
    h.set_ylim(0, min(max(y),1))
    h.title.set_text(key)
    h.set_xlabel('Waiting time in hours')
    h.set_ylabel("Number of occurancies")
    kde = sm.nonparametric.KDEUnivariate(times_dictionary.get(key))
    kde.fit(bw=4, kernel='gau')  # Estimate the densities
    #print(kde.cdf)
    #print(len(kde.density))
    #print(len(kde.support))
    #h.plot(kde.support, kde.density, label="KDE")
    fit_gauss(kde.support, kde.density)
    #y = retrieve_distribution(times_dictionary.get(key))
    #print(y)
    #od = collections.OrderedDict(sorted(y.items()))
    #fit_gauss(list(od.keys()), list(od.values()))
    
    #print(od.values())
    #fit_gauss(od.keys(), od.values())
      
    i += 1
plt.savefig('/Users/a1230101//Documents/GitHub/TimeDistributions/time_plots/DomesticDeclarations.pdf')

"""
Fitting without KDE
"""
"""
fig = plt.figure(figsize=(10,220))
i = 1
for key in sorted(times_dictionary.keys()):
    times = times_dictionary.get(key)
    h = fig.add_subplot(len(times_dictionary), 1, i)
    y_hist, x_hist, _ = h.hist(times, bins=1000, facecolor='g', density = True,)
    h.set_ylim(0, min(max(y_hist),1))
    h.title.set_text(key)
    h.set_xlabel('Waiting time in hours')
    h.set_ylabel("Number of occurancies")
    dict_of_times = retrieve_distribution(times)
    odered_times = collections.OrderedDict(sorted(dict_of_times.items()))
    fit_gauss(x_hist[1:], y_hist)

    #h.plot(kde.support, kde.density, label="KDE")
      
    i += 1
#plt.savefig('/Users/a1230101//Documents/GitHub/TimeDistributions/time_plots/DomesticDeclarationsFit.pdf')
"""