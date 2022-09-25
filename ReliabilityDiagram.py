#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: Naveen Goutham, Camille Le Coz
"""

import numpy as np
import statsmodels.api as sm

class ReliabilityDiagram:
    def __init__(self,observation,forecast,climatology,event_lbound,event_ubound,closed_ends='both',nbins=5,weights=None):
        '''
        Parameters
        ----------
        observation : numpy.ndarray
            Timeseries of observation/truth (1D).
            
        forecast : numpy.ndarray
            Forecast array (2D) with first dimension equal to the dimension of observation,
            and second dimension equal to the ensemble size.
            
        climatology : numpy.ndarray
            Climatological array (2D) with first dimension equal to the dimension of observation, 
            and second dimension equal to the number of years considered in climatology.

        event_lbound : int or float
            The lower bound for the event considered. 
            Example 1: If the event considered is lower tercile, then event_lbound value is 0. 
            Example 2: If the event considered is upper tercile, then event_lbound value is 2/3. 
            NOTE: The value should be in the range of 0 to 1. This value should be lesser than the value of event_ubound.
            
        event_ubound : int or float
            The upper bound for the event considered.
            Example 1: If the event considered is lower tercile, then event_ubound value is 1/3.
            Example 2: If the event considered is upper tercile, then event_ubound value is 1. 
            NOTE: The value should be in the range of 0 to 1. This value should be greater than the value of event_lbound.
            
        closed_ends : str, optional
            The bounds (upper and/or lower) to include in the event formulation. Options: 'left', 'right', 'none', 'both'. 
            Example 1: If closed_ends = 'left', then the event becomes event_lbound <= event < event_ubound. 
            Example-2: If closed_ends = 'both', then the event becomes event_lbound <= event <= event_ubound.
            NOTE: The default value is 'both'. The value is case sensitive.
            
        nbins : int, optional
            Number of bins to stratify the forecasts into. The default is 5.
            NOTE: The bins are of equal width. The number of bins should be lesser than the dimension of observation.
            
        weights : numpy.ndarray, optional
            The weights for the forecast data. It should have the same shape as the forecast data.

        Returns
        -------
        None.

        '''
        #original attributes
        self.ob = observation
        self.fc = forecast
        self.cl = climatology.sort(axis=1)
        self.lb = event_lbound
        self.ub = event_ubound
        self.ends = closed_ends
        self.bins = np.arange(0,1,1/nbins)
        self.weights = weights
        
        #constructed attributes
        self.nsim = len(self.ob)
        self.mem_fc = self.fc.shape[1]
        self.mem_cl = self.cl.shape[1]
            
    def __check_conformity(self):
        '''
        This function checks if every parameter adheres to the basic requirements of the package. This is a private function."
        '''
        if (self.ob.ndim != 1) and (self.fc.ndim and self.cl.ndim != 2) and (self.fc.shape[0] and self.cl.shape[0] != len(self.ob)) and (not (0 <= self.lb and self.ub <= 1)):
            raise ValueError('Please make sure that the input parameters follow the program requirements!')
            exit
        
        if self.lb >= self.ub:
            raise ValueError("Please make sure that event_lbound is lesser than event_ubound.")
            exit
                             
        if len(self.bins) > self.nsim:
            raise ValueError("The number of bins should be lesser than the dimension of observation.")
            exit
      
        list_ends = ["left","right","none","both"]
        if self.ends not in list_ends:
            raise ValueError("Please give a valid entry.")
            exit    
            
        if self.weights is not None and self.weights.shape != self.fc.shape:
            raise ValueError("The shapes of forecast and weights must be the same.")
            exit
        return 
    
    def __get_observed_event(self,l,u):
        '''
        This function checks if the observation has occurred within the given bounds, 
        and returns 1 if the observation has occurred or 0 otherwise.
        '''
        if self.ends == "both":
            return np.logical_and(np.less_equal(l,self.ob), np.less_equal(self.ob,u))
        elif self.ends == "none":
            return np.logical_and(np.less(l,self.ob), np.less(self.ob,u))
        elif self.ends == "left":
            return np.logical_and(np.less_equal(l,self.ob), np.less(self.ob,u))
        elif self.ends == "right":
            return np.logical_and(np.less(l,self.ob), np.less_equal(self.ob,u))
    
    def __get_forecast_probability(self,l,u):
        '''
        This function computes and returns the forecast probability for the defined event by counting the members which are
        within the event definition.
        '''
        if self.weights is None:
            if self.ends == "both":
                return np.sum(((self.fc >= np.tile(l,(self.mem_fc,1)).T)*(self.fc <= np.tile(u,(self.mem_fc,1)).T)), axis=1)/self.mem_fc
            elif self.ends == "none":
                return np.sum(((self.fc > np.tile(l,(self.mem_fc,1)).T)*(self.fc < np.tile(u,(self.mem_fc,1)).T)), axis=1)/self.mem_fc
            elif self.ends == "left":
                return np.sum(((self.fc >= np.tile(l,(self.mem_fc,1)).T)*(self.fc < np.tile(u,(self.mem_fc,1)).T)), axis=1)/self.mem_fc
            elif self.ends == "right":
                return np.sum(((self.fc > np.tile(l,(self.mem_fc,1)).T)*(self.fc <= np.tile(u,(self.mem_fc,1)).T)), axis=1)/self.mem_fc
        else:
            weights = (self.weights.T/np.sum(self.weights, axis=1)).T
            if self.ends == "both":
                return np.sum(((self.fc >= np.tile(l,(self.mem_fc,1)).T)*(self.fc <= np.tile(u,(self.mem_fc,1)).T))*weights, axis=1)/self.mem_fc
            elif self.ends == "none":
                return np.sum(((self.fc > np.tile(l,(self.mem_fc,1)).T)*(self.fc < np.tile(u,(self.mem_fc,1)).T))*weights, axis=1)/self.mem_fc
            elif self.ends == "left":
                return np.sum(((self.fc >= np.tile(l,(self.mem_fc,1)).T)*(self.fc < np.tile(u,(self.mem_fc,1)).T))*weights, axis=1)/self.mem_fc
            elif self.ends == "right":
                return np.sum(((self.fc > np.tile(l,(self.mem_fc,1)).T)*(self.fc <= np.tile(u,(self.mem_fc,1)).T))*weights, axis=1)/self.mem_fc
            
    def contingency_table(self):
        '''
        This function computes and returns the contingency table for the defined event.
        The returned table has "nbins" rows and two columns.  
        
        Returns:
        Contingency table: numpy.ndarray of shape (nbins,2).
            NOTE: The first column corresponds to the "yes" event, whereas the second column corresponds to the "no" event.
        '''
        # Build contingency table
        self.__check_conformity()
        event_names = ['yes', 'no']
        #getting boundaries of the event from climatology
        l, u = np.percentile(self.cl,[100*self.lb,100*self.ub],axis=1)
        if self.lb == 0.0:
            l = -np.inf
        if self.ub == 1.0:
            u = np.inf
        #checking if the event occurred in reality
        event_obs = self.__get_observed_event(l,u)
        #retrieving the probability of occurrence of the event from forecasts
        prob_fc = self.__get_forecast_probability(l,u)
        #getting the bin index for the table
        ind = np.searchsorted(self.bins,prob_fc,side='right') - 1

        #initialise contingency table
        c_table = np.zeros([len(self.bins),len(event_names)])
        #event "yes"
        c_table[:,0] = np.bincount(ind[event_obs],minlength=len(self.bins))
        #event "no"
        c_table[:,1] = np.bincount(ind[~event_obs],minlength=len(self.bins))
        return c_table
    
    def observed_frequency(self):
        '''
        This function returns the observed relative frequency. This is required for plotting the reliability diagram.
        
        Returns:
        Observed relative frequency
        '''
        # Return the observed relative frequency
        self.__check_conformity()
        c_table = self.contingency_table()
        # observed frequency
        return c_table[:,0] / (c_table[:,0] + c_table[:,1])

    def confidence_intervals(self):
        '''
        This function returns the lower and the upper confidence intervals, respectively, corresponding to the observed frequency.
        This is required for plotting the reliability diagram.
        
        Returns:
        lower confidence intervals, upper confidence intervals
        '''
        # Return the confidence intervals
        self.__check_conformity()
        c_table = self.contingency_table()
        # confidence intervals
        ci_low, ci_upp = sm.stats.proportion_confint(c_table[:,0],(c_table[:,0] + c_table[:,1]))
        return ci_low, ci_upp
    
    def forecast_attributes(self):
        '''
        This method computes the brier score, the reliability and the resolution attributes.
        
        NOTE:
        Perfect forecasts will have a brier score of "zero".
        The smaller the reliability component, and the larger the resolution component, the more accurate are the forecasts.
        
        Returns:
        Brier score, Reliability, Resolution
        '''
        self.__check_conformity()
        c_table = self.contingency_table()
        #
        yi = self.bins  # forecast probability
        oi = self.observed_frequency    # observed frequency
        wti = np.sum(c_table,axis=1)/np.sum(c_table)    # weights=number of forecasts yi / total number of forecasts
        om = np.sum(c_table[:,0])/np.sum(c_table)   # overall (unconditional) relative frequency
        o_bar = np.repeat(om,len(self.bins))
        #
        rel = np.nansum(((yi - oi)**2)*wti)    # reliability component of the Brier score
        res = np.nansum(((oi - o_bar)**2)*wti) # resolution component of the Brier score
        bs = rel - res + om*(1-om)             # brier score
        return bs, rel, res
