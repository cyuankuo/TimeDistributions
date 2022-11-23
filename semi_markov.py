import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special
from mult_gauss import MultiGauss
from gauss import Gauss
import math
from convolution import mult_gauss_convolution, mult_gauss_sum, mult_gauss_self_convolution

class SemiMarkov:
    def __init__(self, states, transitions):
        self.states = states 
        # Each transition is a tuple: (from, to, probability, multi Gauss)
        self.transitions = transitions 

    def reduce_node(self, state, label):
        if ((state == 'start') or (state == 'end')):
            return
        else:
            # Calsulate self-loop time
            print("Start calculating self-loop time...")
            for (state, state, _, _) in self.transitions:
                self_loop_time = self.calculate_self_loop_time(state, 0.0001)
            print("End calculating self-loop time...")
            #  Add new transitions
            in_transitions = self.get_in_transitions(state)
            out_transitions = self.get_out_transitions(state)
            number_of_new_transitions = len(in_transitions) * len(out_transitions)
            i = 1
            for in_transition in in_transitions:
                in_state = in_transition[0]
                for out_transition in out_transitions:
                    out_state = out_transition[1]
                    print("Adding transiton " + str(i) + " out of " + str(number_of_new_transitions))
                    i += 1
                    p = self.get_probability(in_state, out_state)
                    time = self.get_time(in_state, out_state)
                    new_p = self.get_probability(in_state,state)*self.get_probability(state,out_state)/(1-self.get_probability(state,state))
                    all_p = p + new_p
                    m1 = self.get_time(in_state, state)
                    m2 = self.get_time(state, out_state)
                    new_time = mult_gauss_convolution(m1,self_loop_time)
                    new_time = mult_gauss_convolution(new_time, m2)
                    all_time = mult_gauss_sum(time, new_time, p/(p+new_p), new_p/(p+new_p))

                    # Remove old transition
                    transition_to_remove = set()
                    for transition in self.transitions:
                        if ((transition[0] == in_state) and (transition[1] == out_state)):
                            transition_to_remove.add(transition)
                    for transition in transition_to_remove:
                        self.transitions.remove(transition) 
                    
                    # Add new transition
                    self.transitions.add((in_state, out_state, all_p, all_time))

                    #all_time.plot_mult_gauss(range(0,1000,1))
                    #print(all_time.probabilities)
                    #plt.title(label=label)
                    #plt.show()

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
        p_current = p
        k = 1
        while (p_current > threshold):
            print(p_current)
            m = mult_gauss_sum(m, mult_gauss_self_convolution(m1, k), 1, p_current)
            p_current *= p
            k += 1
        return m
    
    def get_in_transitions(self, state):
        in_transitions = set()
        for transition in self.transitions:
            if transition[1] == state:
                in_transitions.add(transition)
        return in_transitions
    
    def get_out_transitions(self, state):
        out_transitions = set()
        for transition in self.transitions:
            if transition[0] == state:
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


    