# -*- coding: utf-8 -*-
from assimulo.problem import Implicit_Problem #Imports the problem formulation from Assimulo
import time, math

SEGMENTS = 1

# Extend Assimulos problem definition
class ExtProblem(Implicit_Problem):
    def __init__(self, L_0, V0):
        self.L_0 = L_0
        self.V0  = V0
    # Responsible for handling the events.
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        event_info = event_info[0] # We only look at the state events information.
        while True: #Event Iteration
            self.event_switch(solver, event_info) # Turns the switches            
            b_mode = self.state_events(solver.t, solver.y, solver.yd, solver.sw)
            self.init_mode(solver)     # Pass in the solver to the problem specified init_mode
            a_mode = self.state_events(solver.t, solver.y, solver.yd, solver.sw)            
            event_info = self.check_eIter(b_mode, a_mode)                
            if not True in event_info: # Breaks the iteration loop
                break     
    
    # Helper function for handle_event
    def event_switch(self, solver, event_info):
        """
        Turns the switches.
        """
        length = self.L_0 + (self.V0 * solver.t / SEGMENTS)   
        for i in range(len(event_info)): # Loop across all event functions            
            if event_info[i] != 0:
                pos_ix = 3 * i  + 3      # position index of mass i + 1
                last_pos_ix = pos_ix - 3        
                # calculate the norm of the vector from mass1 to mass0 minus the initial segment length
                solver.sw[i] = math.sqrt((solver.y[pos_ix]-solver.y[last_pos_ix])**2 
                                 + (solver.y[pos_ix+1]-solver.y[last_pos_ix+1])**2 \
                                 + (solver.y[pos_ix+2]-solver.y[last_pos_ix+2])**2) >= length   
        
    # Helper function for handle_event
    def check_eIter(self, before, after):
        """
        Helper function for handle_event to determine if we have event
        iteration.
        
            Input: Values of the event indicator functions (state_events)
            before and after we have changed mode of operations.
        """        
        eIter = [False] * len(before)        
        for i in range(len(before)):
            if (before[i] < 0.0 and after[i] > 0.0) or (before[i] > 0.0 and after[i] < 0.0):
                eIter[i] = True                
        return eIter       
        
    def init_mode(self, solver):
        """
        Initialize the DAE with the new conditions.  """
        solver.make_consistent('IDA_YA_YDP_INIT') # Calculate new initial conditions.
                                                   # see SUNDIALS IDA documentation
                                                   # on the option 'IDA_YA_YDP_INIT'         