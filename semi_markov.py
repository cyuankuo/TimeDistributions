"""
@author: akalenkova (anna.kalenkova@adelaide.edu.au)
"""

from mult_gauss import MultiGauss
from gauss import Gauss
from convolution import mult_gauss_convolution, mult_gauss_sum, threshold

def mean_time_between_events(e1,e2,skip_events,log):
        times = set()
        for trace in log:
            for i in range(len(trace) - 1):
                if trace[i]['concept:name'] == e1:
                    for k in range (i+1, len(trace)):
                        if trace[k]['concept:name'] == e2:
                            time = trace[k]['time:timestamp'] - trace[i]['time:timestamp']
                            times.add(time.total_seconds()//3600)
                            break
                        if not (trace[k]['concept:name'] in skip_events):
                            break
        if len(times) > 0:
            mean = sum(times)/len(times)
        else: 
            mean = 0
        return mean

class SemiMarkov:
    def __init__(self, states, transitions):
        self.states = states 
        # Each transition is a tuple: (from, to, probability, multi Gauss)
        self.transitions = transitions 

    

    def reduce_node(self, state):
        if ((state == 'start') or (state == 'end')):
            return
        else:
            # Calsulate self-loop time
            self_loop_time = MultiGauss([1], [Gauss(0,0)])
            #for transition in self_loops:
            self_loop_time = self.calculate_self_loop_time(state, threshold)
            #  Add new transitions
            in_transitions = self.get_in_transitions(state)
            out_transitions = self.get_out_transitions(state)
            i = 1
            for in_transition in in_transitions:
                in_state = in_transition[0]
                for out_transition in out_transitions:
                    out_state = out_transition[1]
                    p = self.get_probability(in_state, out_state)
                    time = self.get_time(in_state, out_state)
                    new_p = self.get_probability(in_state,state)*self.get_probability(state,out_state)/(1-self.get_probability(state,state))
                    all_p = p + new_p
                    m1 = self.get_time(in_state, state)
                    m2 = self.get_time(state, out_state)
                    new_time = mult_gauss_convolution(m1, self_loop_time)
                    new_time = mult_gauss_convolution(new_time, m2)
                    all_time = mult_gauss_sum(time, new_time, p/all_p, new_p/all_p)

                    # Remove old transition
                    transition_to_remove = set()
                    for transition in self.transitions:
                        if ((transition[0] == in_state) and (transition[1] == out_state)):
                            transition_to_remove.add(transition)
                    for transition in transition_to_remove:
                        self.transitions.remove(transition)
                            
                    # Add new transition
                    self.transitions.add((in_state, out_state, all_p, all_time))


            # Remove state
            transition_to_remove = set()
            for transition in self.transitions:
                if ((transition[0] == state) or (transition[1] == state)):
                    transition_to_remove.add(transition)
            for transition in transition_to_remove:
                 self.transitions.remove(transition)
            self.states.remove(state)


    def  calculate_self_loop_time(self, state, threshold):
        m1 = self.get_time(state, state)
        p = self.get_probability(state, state)
        m = MultiGauss([1-p],[Gauss(0,0)])
        p_current = p * (1-p)
        conv = MultiGauss([1],[Gauss(0,0)])
        while (p_current > threshold):
            conv = mult_gauss_convolution(m1, conv)
            m = mult_gauss_sum(m, conv, 1, p_current)
            p_current *= p
        m.normalise_gauss()
        return m
    
    def get_in_transitions(self, state):
        in_transitions = set()
        for transition in self.transitions:
            if (transition[1] == state) and (transition[0] != state):
                in_transitions.add(transition)
        return in_transitions
    
    def get_in_transitions_with_loop(self, state):
        in_transitions = set()
        for transition in self.transitions:
            if (transition[1] == state):
                in_transitions.add(transition)
        return in_transitions
    
    def get_out_transitions(self, state):
        out_transitions = set()
        for transition in self.transitions:
            if (transition[0] == state) and (transition[1] != state):
                out_transitions.add(transition)
        return out_transitions
    
    def get_out_transitions_with_loop(self, state):
        out_transitions = set()
        for transition in self.transitions:
            if (transition[0] == state):
                out_transitions.add(transition)
        return out_transitions
 
    def get_probability(self, state1, state2):
        for transition in self.transitions:
            if ((transition[0] == state1) and (transition[1] == state2)):
                return transition[2]
        return 0
    
    def get_time(self, state1, state2):
        for transition in self.transitions:
            if ((transition[0] == state1) and (transition[1] == state2)):
                return transition[3]
        return MultiGauss([],[])
    
    def select_next_state(self):
        min_degree = 100000000
        for state in self.states:
            degree = 0
            degree = len(self.get_in_transitions_with_loop(state))*len(self.get_out_transitions_with_loop(state))
            if (degree < min_degree) and (state != 'start') and (state != 'end'):
                next_state = state
                min_degree = degree
        return next_state



    