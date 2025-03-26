"""
@author: akalenkova (anna.kalenkova@adelaide.edu.au)
"""

from pm4py.objects.log.importer.xes import importer as xes_importer
import matplotlib
matplotlib.use('Agg')
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
import numpy as np, scipy.stats as st
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
from datetime import datetime
from pm4py.objects.log.obj import EventLog

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

def extract_times_event_log(log):
    times = []
    print(f"log length is {len(log)}")
    for trace in log:
        start = trace[0]['time:timestamp']
        end = trace[len(trace)-1]['time:timestamp']
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

    # print(f"states is: {states}")
    # print(f"dfg is: {dfg}")

    for key1 in states:
        out_frequences[key1] = 0
        for key2 in states:
            if (key1, key2) in dfg and dfg[key1,key2] > 0:
                out_frequences[key1] += dfg[key1,key2]

    # print(f"multi_gausses is: {multi_gausses.keys()}")
    for key1 in states:
        for key2 in states:
            if (key1, key2) in dfg and dfg[key1,key2] > 0:

                if ((key1 == 'start') or  (key1 == 'end') or (key2 == 'start') or  (key2 == 'end')):
                    transitions.add(tuple([key1, key2, dfg[key1,key2]/out_frequences[key1], MultiGauss([1], [Gauss(0, 0)])]))
                else:
                    transitions.add(tuple([key1, key2, dfg[key1,key2]/out_frequences[key1],
                    multi_gausses["['" + str(key1) + "', '" + str(key2) + "']"]]))
    print()
    print('DFG is built')
    return SemiMarkov(states, transitions)

def process_transition(args, bw=4, counter=None):
    key, times = args
    try:
        kde = sm.nonparametric.KDEUnivariate(times)
        kde.fit(bw=bw, kernel='gau')
        multi_gauss = fit_gauss(kde.support, kde.density, times)
        if counter is not None:
            counter.value += 1
        print(f"{counter.value} transitions fitted")
        return str([key.partition('->')[0], key.partition('->')[2]]), multi_gauss
    except Exception as e:
        print(f"Error processing {key}: {str(e)}")
        return str([key.partition('->')[0], key.partition('->')[2]]), MultiGauss([1], [Gauss(0, 0)])

if __name__ == '__main__':
    # move print statements to log
    log_file = open('filtered_hospital.log', 'a')
    sys.stdout = log_file

    print("================Start of Log================")



    variant = xes_importer.Variants.ITERPARSE
    parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
    log = xes_importer.apply('logs/' + sys.argv[1], variant=variant, parameters=parameters)
    discovery_times = {}
    fitting_times = {}
    reduction_times = {}
    kl_divergences = {}
    filtered_traces = EventLog()
    # filter log for 2FC
    for i, trace in enumerate(log):
        if any(event.get('concept:name') == "2FC" for event in trace):
            filtered_traces.append(trace)

    print(f"prev total activities: {len(log)}")
    print(f"total 2FC activities: {len(filtered_traces)}")
    event_log_times = extract_times_event_log(filtered_traces)

    filtered_event_log_times = []
    for times in event_log_times:
        if times < 1200:
            filtered_event_log_times.append(times)

    print(f" filtered event log times is :{len(filtered_event_log_times)}")
    for k in [1,2]:

        print()
        print("Order: k=" + str(k))

        start = time.time()
        # log_for_discovery = deepcopy(log)
        log_for_discovery = deepcopy(filtered_traces)
        times_dictionary = {}
        log_processed = log_parser.prepare_log(log_for_discovery, k)
        print(log_processed)
        dfg, start_activities, end_activities = discover_dfg(log_processed)
        dfg["end", "start"] = 1

        filtered_dfg = {}
        for (source, target), freq in dfg.items():
            source_contains_2FC = source == "2FC" or "2FC," in source or ",2FC" in source or ",2FC," in source
            target_contains_2FC = target == "2FC" or "2FC," in target or ",2FC" in target or ",2FC," in target

            if (source_contains_2FC or target_contains_2FC or
                    (source == 'start' and target == '2FC') or  # Only start -> 2FC
                    (target == 'end' and source == '2FC')):  # Only 2FC -> end
                filtered_dfg[(source, target)] = freq
        # Ensure cycle back if needed
        if any('2FC' in s for s in filtered_dfg.keys()):  # Only add if 2FC exists
            filtered_dfg[('end', 'start')] = 1

        print(f"filtered dfg is: {filtered_dfg}")
        end = time.time()
        print()
        print("Discovery time:")
        print(end-start)
        if k not in discovery_times:
            discovery_times[k] = {end-start}
        else:
            discovery_times[k].add(end-start)

        "Express analysis"
        # cut the log to get better precision of the limiting probabilities
        number_of_chunks = len(log_for_discovery)
        overall_times = []
        #cnt = 0
        temp_log_for_discovery = deepcopy(log_for_discovery)
        for traces in np.array_split(np.array(temp_log_for_discovery, dtype=object), number_of_chunks):
            processed_traces = log_parser.prepare_log(traces, k)
            #for i in range(len(traces[0])):
            #    print(traces[0][i])
            #print(cnt)
        #    cnt += 1
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

        print(str(round(estimated_mean_time//86400)) + 'd ' + str(round(estimated_mean_time%86400//3600)) + 'h ' + str(round(estimated_mean_time%3600//60)) + 'm ' + str(round(estimated_mean_time%60)) + 's ')


        "Full analysis"

        start = time.time()
        extract_times_with_future(log_processed)

        """
        Fitting using Gaussian KDE
        """
        mult_gausses = {}
        filtered_times_dictionary = {key: value for key, value in times_dictionary.items() if "2FC->" in key or "->2FC" in key or ",2FC" in key or "2FC," in key}
        # total_transitions = len(times_dictionary)
        total_filtered_transitions = len(filtered_times_dictionary)

        # print(times_dictionary.keys())
        print(f"filtered dict keys is: {filtered_times_dictionary.keys()}")
        # transition_data = [(key, times_dictionary[key]) for key in sorted(times_dictionary.keys())]
        transition_data = [(key, times_dictionary[key]) for key in sorted(filtered_times_dictionary.keys())]

        n_processes = max(1, cpu_count() - 1)
        print(f"Using {n_processes} processes", file=log_file)
        log_file.flush()

        # Create a shared counter
        manager = Manager()
        counter = manager.Value('i', 0)

        # print(f"Total transitions to process: {total_transitions}", file=log_file)
        print(f"Total transitions to process: {total_filtered_transitions}", file=log_file)

        log_file.flush()

        with Pool(processes=n_processes) as pool:
            print("Pool created, mapping transitions...", file=log_file)
            log_file.flush()
            process_func = partial(process_transition, bw=4, counter=counter)
            results = pool.map(process_func, transition_data)


            for i, (key, multi_gauss) in enumerate(results):
                mult_gausses[key] = multi_gauss
                # if (i + 1) % max(1, total_transitions // 10) == 0:
                #     print(f"{counter.value} out of {total_transitions} transitions fitted (iteration {i + 1})", file=log_file)
                #     log_file.flush()
                if (i + 1) % max(1, total_filtered_transitions // 10) == 0:
                    print(f"{counter.value} out of {total_filtered_transitions} transitions fitted (iteration {i + 1})", file=log_file)
                    log_file.flush()



        # for key in sorted(times_dictionary.keys()):
        #     # print(f"Input data: {times_dictionary.get(key)}")
        #     kde = sm.nonparametric.KDEUnivariate(times_dictionary.get(key))
        #     kde.fit(bw=4, kernel='gau')  # Estimate the densities
        #     # print(f"Support min: {kde.support.min()}, max: {kde.support.max()}")
        #     # print(f"Density min: {kde.density.min()}, max: {kde.density.max()}")
        #     # print(f"Support: {kde.support[:5]}... (length: {len(kde.support)})")
        #     # print(f"Density: {kde.density[:5]}... (max: {kde.density.max()})")
        #     multi_gauss = fit_gauss(kde.support, kde.density, times_dictionary.get(key))
        #     mult_gausses[str([key.partition('->')[0], key.partition('->')[2]])] = multi_gauss
        #     # count += 1
        #     # print(f"{count} out of {total_transitions} transitions fitted.......")

        end = time.time()
        print("Fitting time:", file=log_file)
        print(end - start, file=log_file)
        log_file.flush()
        if k not in fitting_times:
            fitting_times[k] = {end-start}
        else:
            fitting_times[k].add(end-start)

        # semi_markov = build_semi_markov(dfg, mult_gausses)
        print(f"filtered dfg is: {filtered_dfg}")
        semi_markov = build_semi_markov(filtered_dfg, mult_gausses)



        print("Number of states: " + str(len(semi_markov.states)))
        print("Number of transitions: " + str(len(semi_markov.transitions)))
        state_degrees = semi_markov.state_degrees()
        avg_state_degree = np.average(state_degrees)
        print("Average state degree: " + str(avg_state_degree))
        max_state_degree = np.max(state_degrees)
        print("Max state degree: " + str(max_state_degree))


        print("Simulation...")
        times = semi_markov.simulate()
        print(times)

        states = deepcopy(semi_markov.states)
        print(f"semi_markov states: {len(states)}")
        start = time.time()
        while len(semi_markov.states) > 2:
            next_state = semi_markov.select_next_state()
            semi_markov.reduce_node(next_state)
        end = time.time()
        print()
        print("Reduction time:")
        print(end-start)
        if k not in reduction_times:
            reduction_times[k] = {end-start}
        else:
            reduction_times[k].add(end-start)

        for idx, transition in enumerate(semi_markov.transitions):
            print(f"current transition is: {transition}")
            if transition[0] == 'start':
                multi_gauss = transition[3]
                multi_gauss.remove_zero()
                color = {
                    1: "tab:red",
                    2: "tab:blue",
                    3: "k",
                    4: "tab:green",
                    5: "tab:purple"
                }
                multi_gauss.plot_trunc_mult_gauss(range(-10,400,1), label="Semi-Markov Model, order="+str(k), color = color.get(k))


                print()
                print("Peaks:")
                print(multi_gauss.calculate_peaks())

                print()
                print("KL Divergence:")
                kl_divergence = multi_gauss.calc_kl_divergence(20, filtered_event_log_times, event_log_times)
                print(kl_divergence)
                if k not in kl_divergences:
                    kl_divergences[k] = {kl_divergence}
                else:
                    kl_divergences[k].add(kl_divergence)
                print()


    for k in [1,2]:
        print()
        print("Metrics for order " + str(k) + ":")
        discovery_times_values = discovery_times[k]
        print(discovery_times_values)
        discovery_times_average = np.mean(list(discovery_times_values))
        print("Discovery times average:")
        print(discovery_times_average)
        discovery_times_interval = st.t.interval(0.95, df=len(discovery_times_values)-1,
                  loc=np.mean(list(discovery_times_values)),
                  scale=st.sem(list(discovery_times_values)))
        print("Discovery times interval:")
        print(discovery_times_interval[1]-discovery_times_average)
        print()

        fitting_times_values = fitting_times[k]
        print(fitting_times_values)
        fitting_times_average = np.mean(list(fitting_times_values))
        print("Fitting times average:")
        print(fitting_times_average)
        fitting_times_interval = st.t.interval(0.95, df=len(fitting_times_values)-1,
                  loc=np.mean(list(fitting_times_values)),
                  scale=st.sem(list(fitting_times_values)))
        print("Fitting times interval:")
        print(fitting_times_interval[1]-fitting_times_average)
        print()

        reduction_times_values = reduction_times[k]
        print(reduction_times_values)
        reduction_times_average = np.mean(list(reduction_times_values))
        print("Reduction times average:")
        print(reduction_times_average)
        reduction_times_interval = st.t.interval(0.95, df=len(reduction_times_values)-1,
                  loc=np.mean(list(reduction_times_values)),
                  scale=st.sem(list(reduction_times_values)))
        print("Reduction times interval:")
        print(reduction_times_interval[1]-reduction_times_average)
        print()

        kl_divergence_values = kl_divergences[k]
        print(kl_divergence_values)
        kl_divergence_average = np.mean(list(kl_divergence_values))
        print("KL-divergence average:")
        print(kl_divergence_average)
        kl_divergence_interval = st.t.interval(0.95, df=len(kl_divergence_values)-1,
                  loc=np.mean(list(kl_divergence_values)),
                  scale=st.sem(list(kl_divergence_values)))
        print("KL-divergence interval:")
        print(kl_divergence_interval[1]-kl_divergence_average)
        print()

        sys.stdout = sys.__stdout__
        log_file.close()


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
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    plt.savefig(f'event_log_histogram_{timestamp}.png')
    # plt.savefig('event_log_histogram.png')
    plt.show()


