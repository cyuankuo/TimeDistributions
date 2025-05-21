import pm4py
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery
from pm4py.visualization.dfg import visualizer as dfg_visualizer
from pm4py.algo.filtering.dfg import dfg_filtering
from pm4py.statistics.attributes.log.get import get_attribute_values
from pm4py.objects.log.importer.xes import importer as xes_importer
from datetime import datetime, timezone
from collections import defaultdict


# Load the event log
file_path = ""
# event_log = pm4py.read_xes(file_path)
event_log = xes_importer.apply(file_path)

def build_intervals_from_event_log(event_log):
    ward_intervals = defaultdict(list)

    # print(list(event_log[0][0].keys()))

    for trace in event_log:
        patient = ""
        events = []
        for event in trace:
            timestamp = event.get("time:timestamp", None)
            ward = event.get("ward", None)
            patient = event.get("URN") if patient == "" else patient
            if ward is not None and timestamp is not None:
                events.append({
                    "timestamp" : timestamp,
                    "ward": ward})
        for i, event in enumerate(events):
            ward = event["ward"]
            enter_time = event["timestamp"]
            exit_time = events[i + 1]["timestamp"] if i + 1 < len(events) else None
            ward_intervals[ward].append({
                "patient_id": patient,
                "enter": enter_time,
                "exit": exit_time
            })
    return ward_intervals


# Function to get occupancy from the event log
def current_occupancy(timestamp, ward, ward_intervals):
    count = 0
    for interval in ward_intervals.get(ward, []):
        if interval["enter"] <= timestamp and (interval["exit"] is None or timestamp < interval["exit"]):
            count += 1
    return count

ward_intervals = build_intervals_from_event_log(event_log)

timestamp = ts = datetime(2011, 1, 11, 11, 11, 0, tzinfo=timezone.utc)
ward = ""
print(f"current_occupancy for ward {ward} at {timestamp.strftime('%Y-%m-%d %H:%M:%S')}: {current_occupancy(timestamp, ward, ward_intervals)}")


# wards_with_multiple_entries = {
#     ward: intervals
#     for ward, intervals in ward_intervals.items()
#     if len(intervals) > 1
# }
#
# min_ward = None
# min_count = None
#
# for ward, intervals in wards_with_multiple_entries.items():
#     entry_count = len(intervals)
#
#     if min_count is None or entry_count < min_count:
#         min_count = entry_count
#         min_ward = ward
#
#
# print(f"Ward with >1 entry and minimum count: {min_ward} ({min_count} entries)")


#################
#################
# # Discover the DFG
# result = dfg_discovery.apply(event_log, variant=dfg_discovery.Variants.FREQUENCY)
# print(f"DFG result: {result}")
# if isinstance(result, tuple):
#     dfg, start_activities, end_activities = result
# else:
#     dfg = result
#     start_activities = pm4py.get_start_activities(event_log)
#     end_activities = pm4py.get_end_activities(event_log)
#
# # Extract activity counts
# activities_count = get_attribute_values(event_log, "concept:name")
#
# # Debug: Print input data
# print(f"Original DFG: {dfg}")
# print(f"Start Activities: {start_activities}")
# print(f"End Activities: {end_activities}")
# print(f"Activities Count: {activities_count}")
#
# # Filter DFG with activities (top 50%)
# filtered_result = dfg_filtering.filter_dfg_on_activities_percentage(
#     dfg, start_activities, end_activities, activities_count, percentage=0.2
# )
# print(f"Filtered result type: {type(filtered_result)}")
# print(f"Filtered result length: {len(filtered_result)}")
# print(f"Filtered result content: {filtered_result}")
#
# # Extract components from the tuple
# filtered_dfg, filtered_start_activities, filtered_end_activities, filtered_activities_count = filtered_result
# print(f"Extracted DFG: {filtered_dfg}")
# print(f"Filtered Start Activities: {filtered_start_activities}")
# print(f"Filtered End Activities: {filtered_end_activities}")
# print(f"Filtered Activities Count: {filtered_activities_count}")
#
# # Check if DFG is empty
# if not filtered_dfg:
#     raise ValueError("Filtered DFG is empty.")
#
# # Visualize the filtered DFG
# parameters = dfg_visualizer.Variants.FREQUENCY.value.Parameters
# gviz = dfg_visualizer.apply(filtered_dfg, log=event_log, variant=dfg_visualizer.Variants.FREQUENCY,
#                            parameters={
#                                parameters.FORMAT: "png",
#                                parameters.RANKDIR: "TB",
#                                "graph_attr": {"splines": "ortho", "nodesep": "0.5", "ranksep": "1"}
#                            })
#
# # Save and view the visualization
# dfg_visualizer.save(gviz, "dfg_visualization.png")
# dfg_visualizer.view(gviz)