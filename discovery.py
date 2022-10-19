from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.objects.dfg.exporter import exporter as dfg_exporter
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery
from pm4py.statistics.start_activities.log import get as start_activities
from pm4py.statistics.end_activities.log import get as end_activities
from pm4py.visualization.dfg import visualizer as dfg_visualization
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
                if not event['concept:name'] + '_' + next_event['concept:name'] in times_dictionary.keys():
                    times_dictionary[event['concept:name'] + '_' + next_event['concept:name']] = [time.total_seconds()//3600]
                else:
                    times_dictionary[event['concept:name'] + '_' + next_event['concept:name']].append(time.total_seconds()//3600)
            event = next_event
            first = False

variant = xes_importer.Variants.ITERPARSE
parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
log = xes_importer.apply('/Users/a1230101//Documents/GitHub/TimeDistributions/bpi_challenge_2013_incidents.xes', 
    variant=variant, parameters=parameters)
times_dictionary = {}

extract_times()

#fig = plt.figure(figsize=(10,300))
#i = 1
#for key in sorted(times_dictionary.keys()):
    #h = fig.add_subplot(len(times_dictionary), 1, i)
    #print(h)
    #h.hist(times_dictionary.get(key), bins=1000, facecolor='g')
    #h.title.set_text(key)
    #h.set_xlabel('Number of hours')
    #i += 1
plt.hist(times_dictionary.get("Completed"), bins=100, facecolor='g')
plt.xlabel("Waiting time in hours")
plt.ylabel("Number of occurancies")
plt.show()
plt.savefig('completed.eps')

#h.show()

#sa = start_activities.get_start_activities(log)
#ea = end_activities.get_end_activities(log)
#print(sa)
#print(ea)


#dfg = dfg_discovery.apply(log, variant=dfg_discovery.Variants.FREQUENCY)

#for start in sa:
#    for end in ea:
#        dfg[end, start] += sa[start]

#dfg_exporter.apply(dfg, '/Users/a1230101//Documents/GitHub/MCDiscovery/output.dfg', parameters={dfg_exporter.Variants.CLASSIC.value.Parameters.START_ACTIVITIES: sa,
#                                   dfg_exporter.Variants.CLASSIC.value.Parameters.END_ACTIVITIES: ea})

#gviz = dfg_visualization.apply(dfg, log=log, variant=dfg_visualization.Variants.FREQUENCY)
#dfg_visualization.view(gviz)
