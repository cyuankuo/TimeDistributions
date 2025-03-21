"""
@author: akalenkova (anna.kalenkova@adelaide.edu.au)
"""

import numpy as np

from mult_gauss import MultiGauss
from gauss import Gauss
from convolution import mult_gauss_convolution, mult_gauss_sum


self_loop_threshold = 0.1

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

    def start_state(self):
        for state in self.states:
            if state == 'start':
                return state
    
    def end_state(self):
        for state in self.states:
            if state == 'end':
                return state
    

    def draw_transition(self, state):
        transitions = self.get_out_transitions_with_loop(state)
        prob = []
        trans = []

        for transition in transitions:
            prob.append(transition[2])
            trans.append(transition)
        print(prob)
        index = np.random.choice(len(trans), 1, p=prob)
        transition = trans[index[0]]
        return transition

    def draw_time(self, transition):
        multi_gauss = transition[3]
        gauss = np.random.choice(multi_gauss.gaussians, 1, p=multi_gauss.probabilities)
        time = np.random.normal(loc=gauss[0].mean, scale=gauss[0].deviation)
        return time


    def simulate(self):
        times = set()
        iterations = 1
        for i in range(iterations):
            time = 0
            state = self.start_state()
            end = self.end_state() 
            while state != end:
                transition = self.draw_transition(state)
                time += self.draw_time(transition)
                state = transition[1]
                print(state)
            times.add(time)
        return times

    def reduce_node(self, state):
        if ((state == 'start') or (state == 'end')):
            return
        else:
            # Calsulate self-loop time
            self_loop_time = MultiGauss([1], [Gauss(0,0)])
            #for transition in self_loops:
            self_loop_time = self.calculate_self_loop_time(state)
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


    def  calculate_self_loop_time(self, state):
        m1 = self.get_time(state, state)
        p = self.get_probability(state, state)
        m = MultiGauss([1-p],[Gauss(0,0)])
        p_current = p * (1-p)
        conv = MultiGauss([1],[Gauss(0,0)])
        while (p_current > self_loop_threshold):
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

    def state_degrees(self):
        state_degrees = []
        for state in self.states:
            if ((state != 'start') and (state != 'end')):
                in_transitions = len(self.get_in_transitions(state))
                out_transitions = len(self.get_out_transitions(state))
                state_degrees.append(max(in_transitions, out_transitions))
        return state_degrees

    