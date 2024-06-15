"""
@author: akalenkova (anna.kalenkova@adelaide.edu.au)
"""

from pm4py.objects.log.importer.xes import importer as xes_importer
import matplotlib.pyplot as plt
import statsmodels.api as sm
from fit_distribution import fit_gauss
import log_parser
from pm4py import discover_dfg
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery
import time 
from copy import deepcopy
from mult_gauss import MultiGauss
from gauss import Gauss
from semi_markov import SemiMarkov
import sys, dfg_utils, stat_utils
import numpy as np


def extract_times_with_future(log):
    for trace in log:
        first = True
        for next_event in trace:
            if not first and next_event['concept:name'] != 'end' and event['concept:name'] != 'start':
                time = next_event['time:timestamp'] - event['time:timestamp']
                if not event['concept:name'] + '->' + next_event['concept:name'] in times_dictionary.keys():
                    times_dictionary[event['concept:name'] + '->' + next_event['concept:name']] = [time.total_seconds()//3600]
                else:
                    times_dictionary[event['concept:name'] + '->' + next_event['concept:name']].append(time.total_seconds()//3600)
            event = next_event
            first = False

def extract_times_event_log():
    times = []
    for trace in log:
        start = trace[0]['time:timestamp']
        end = trace[len(trace)-1] ['time:timestamp']   
        time = end - start
        times.append(time.total_seconds()//3600)
    return times

# retrieve distribution of values from given y
def retrieve_distribution(y):
    result = {}
    for i in y:
        count_i = y.count(i)
        result[i] = count_i
        result[i] /= len(y)
    return result

def build_semi_markov(dfg, multi_gausses):
    
    states = set()
    transitions = set()
    out_frequences = {}

    for key in dfg.keys():
        states.add(key[0])
    
    for key1 in states:
        out_frequences[key1] = 0
        for key2 in states:
            if dfg[key1,key2] > 0:
                out_frequences[key1] += dfg[key1,key2]


    for key1 in states:
        for key2 in states:
            if dfg[key1,key2] > 0:
                if ((key1 == 'start') or  (key1 == 'end') or (key2 == 'start') or  (key2 == 'end')):
                    transitions.add(tuple([key1, key2, dfg[key1,key2]/out_frequences[key1], MultiGauss([1], [Gauss(0, 0)])]))
                else:
                    transitions.add(tuple([key1, key2, dfg[key1,key2]/out_frequences[key1], 
                    multi_gausses["['" + str(key1) + "', '" + str(key2) + "']"]]))
    print()
    print('DFG is built')
    return SemiMarkov(states, transitions)


variant = xes_importer.Variants.ITERPARSE
parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
log = xes_importer.apply('logs/' + sys.argv[1], variant=variant, parameters=parameters)

event_log_times = extract_times_event_log()

filtered_event_log_times = []
for times in event_log_times:
    if times < 1000:
        filtered_event_log_times.append(times)


for k in [1,2,3,4,5]:
 
    print()
    print("Order: k=" + str(k))

    start = time.time()
    log_for_discovery = deepcopy(log)
    times_dictionary = {}
    log_processed = log_parser.prepare_log(log_for_discovery, k)
    dfg, start_activities, end_activities = discover_dfg(log_processed)
    dfg["end", "start"] = 1

    end = time.time()
    print()
    print("Discovery time:")
    print(end-start)

    "Express analysis"

    # cut the log to get better precision of the limiting probabilities
    number_of_chunks = len(log_for_discovery)
    overall_times = []
    temp_log_for_discovery = deepcopy(log_for_discovery)
    for traces in np.array_split(np.array(log_for_discovery, dtype=object), number_of_chunks):

        dfg_express = dfg_discovery.apply(traces, variant=dfg_discovery.Variants.FREQUENCY)
        dfg_express["end", "start"] = 1
        log_activities=log_parser.log_activities(traces)
        times = log_parser.calculate_times(traces) 
        means = stat_utils.calculate_means(dfg_express, times, log_activities)
        #print(means)
        limiting_probabilities = dfg_utils.calculate_limiting_probabilities(dfg_express, log_activities)
        #print(limiting_probabilities)
        overall_time = 0
        for i in range(0, len(log_activities)):
            overall_time += limiting_probabilities[log_activities[i]]*means[log_activities[i]]
        overall_time /= limiting_probabilities['start']
        overall_times.append(overall_time)
    estimated_mean_time = np.average(overall_times)
  
    print()
    print("Mean time of the process:")
    print(str(round(estimated_mean_time//86400)) + 'd ' + str(round(estimated_mean_time%86400//3600)) + 'h ' + str(round(estimated_mean_time%3600//60)) + 'm ' + str(round(estimated_mean_time%60)) + 's ')

    
    "Full analysis"
    
    start = time.time()
    extract_times_with_future(log_processed)

    """
    Fitting using Gaussian KDE
    """
    mult_gausses = {}
    for key in sorted(times_dictionary.keys()):
        kde = sm.nonparametric.KDEUnivariate(times_dictionary.get(key))
        kde.fit(bw=4, kernel='gau')  # Estimate the densities
        multi_gauss = fit_gauss(kde.support, kde.density, times_dictionary.get(key))
        mult_gausses[str([key.partition('->')[0], key.partition('->')[2]])] = multi_gauss

    end = time.time()
    print()
    print("Fitting time:")
    print(end - start)

    semi_markov = build_semi_markov(dfg, mult_gausses)
    print("Number of states: " + str(len(semi_markov.states)))
    print("Number of transitions: " + str(len(semi_markov.transitions)))
    states = deepcopy(semi_markov.states)
    
    start = time.time()
    while len(semi_markov.states) > 2:
        next_state = semi_markov.select_next_state()
        semi_markov.reduce_node(next_state)
    end = time.time()
    print()
    print("Reduction time:")
    print(end-start)

    for transition in semi_markov.transitions:
        if transition[0] == 'start':
            multi_gauss = transition[3]
            multi_gauss.remove_zero()
            color = {
                1: "tab:red",
                2: "tab:blue",
                3: "k",
                4: "tab:green",
                5: "tab:purple",
                10: "tab:red"
            }
            multi_gauss.plot_trunc_mult_gauss(range(-10,400,1), label="Semi-Markov Model, order="+str(k), color = color.get(k))
            print()
            print("Peaks:")
            print(multi_gauss.calculate_peaks())
                
            print()
            print("KL Divergence:")
            print(multi_gauss.calc_kl_divergence(20, filtered_event_log_times))
            print()
                


"""
Plotting event log
"""

cm = plt.cm.get_cmap('OrRd')
y, x, _ = plt.hist(filtered_event_log_times, bins=150, fc=cm(0.25), density=True, label='Event log')
 
plt.xlim([-10, 400])
plt.legend(loc="upper right")
plt.title('')
plt.xlabel('Overall time in hours')
plt.ylabel('Probability')
plt.show()


