"""
@author: akalenkova (anna.kalenkova@adelaide.edu.au)
"""

import statistics

def calculate_means(dfg,times,log_activities):
    all_means = {}
    for activity1 in log_activities:
        i = 0
        mean = 0
        for activity2 in log_activities:
            if dfg[activity1,activity2] > 0:
                #print(times[activity1+";"+activity2])
                for time in times[activity1+";"+activity2]:
                    mean += time.total_seconds()
                    i += 1
        mean /= i
        all_means[activity1]=mean
    return all_means

def calculate_standard_deviation_times(dfg,times,log_activities):
    all_standard_deviation_times = {}
    all_times = {}
    for activity1 in log_activities:
        for activity2 in log_activities:
            if dfg[activity1,activity2] > 0:
                all_times[activity1] = []
                for time in times[activity1+";"+activity2]:
                    all_times[activity1].append(time.total_seconds())

    for activity1 in all_times:
            times = all_times[activity1]
            if (len(times) > 1):
                all_standard_deviation_times[activity1] = statistics.stdev(times)
            else:
                all_standard_deviation_times[activity1] = 0

    return all_standard_deviation_times

