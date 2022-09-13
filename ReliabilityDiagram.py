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
            NOTE: The value should be in the range of 0 to 1.
            
        event_ubound : int or float
            The upper bound for the event considered.
            Example 1: If the event considered is lower tercile, then event_ubound value is 1/3.
            Example 2: If the event considered is upper tercile, then event_ubound value is 1. 
            NOTE: The value should be in the range of 0 to 1.
            
        closed_ends : str, optional
            The bounds (upper and/or lower) to include in the event formulation. Options: 'left', 'right', 'none', 'both'. 
            Example 1: If closed_ends = 'left', then the event becomes event_lbound <= event < event_ubound. 
            Example-2: If closed_ends = 'both', then the event becomes event_lbound <= event <= event_ubound.
            NOTE: The default value is 'both'. The value is case sensitive.
            
        nbins : int, optional
            Number of bins to stratify the forecasts into. The default is 5.
            NOTE: The bins are of equal width.
            
        weights : numpy.ndarray, optional
            The weights for the forecast data. It should have the same shape as the forecast data.

        Returns
        -------
        None.

        '''
        #original attributes
        self.ob = observation
        self.fc = forecast
        self.cl = np.sort(climatology,axis=1)
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
        if (self.ob.ndim != 1) and (self.fc.ndim and self.cl.ndim != 2) and (self.fc.shape[0] and self.cl.shape[0] != len(self.ob)) and (not (0 <= self.lb and self.ub <= 1)):
            raise ValueError('Please make sure that the input parameters follow the program requirements!')
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
        #add condition for boundaries. i.e., if lbound = 0, then ob <= ubound, and if ubound = 1, then ob >= lbound
        if self.ends == "both":
            return np.logical_and(np.less_equal(l,self.ob), np.less_equal(self.ob,u))
        elif self.ends == "none":
            return np.logical_and(np.less(l,self.ob), np.less(self.ob,u))
        elif self.ends == "left":
            return np.logical_and(np.less_equal(l,self.ob), np.less(self.ob,u))
        elif self.ends == "right":
            return np.logical_and(np.less(l,self.ob), np.less_equal(self.ob,u))
    
    def __get_forecast_probability(self,l,u):
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
    
    def observed_frequency_confidence(self):
        # Return the observed relative frequency and the corresponding confidence intervals
        self.__check_conformity()
        c_table = self.contingency_table()
        #
        self.oi = c_table[:,0] / (c_table[:,0] + c_table[:,1])  # observed frequency
        self.ci_low, self.ci_upp = sm.stats.proportion_confint(c_table[:,0],(c_table[:,0] + c_table[:,1])) # confidence interval
        return self.oi, self.ci_low, self.ci_upp
        
    def forecast_attributes(self):
        # Compute the reliability and resolution components of the Brier score from the contingency table
        self.__check_conformity()
        c_table = self.contingency_table()
        #
        yi = self.bins  # forecast probability
        oi = c_table[:,0] / (c_table[:,0] + c_table[:,1])    # observed frequency
        wti = np.sum(c_table,axis=1)/np.sum(c_table)    # weights=number of forecasts yi / total number of forecasts
        rel = np.nanmean(((yi - oi)**2)*wti)    # reliability component of the Brier score
        om = np.sum(c_table[:,0])/np.sum(c_table)   # overall (unconditional) relative frequency
        o_bar = np.repeat(om,len(self.bins))
        res = np.nanmean(((oi - o_bar)**2)*wti) # resolution component of the Brier score
        #print('Reliability = ',round(rel,4),' | Resolution = ',round(res,4),'\n(Reliability - Resolution) = ',round(rel - res,4))
        return rel, res
        
#    def plot_diagram(self):
#        # Plot reliability diagram
#        self.__check_conformity()
#        c_table = self.contingency_table()
#        #
#        p = c_table[:,0] / (c_table[:,0] + c_table[:,1])    # observed frequency
#        ci_low, ci_upp = sm.stats.proportion_confint(c_table[:,0],(c_table[:,0] + c_table[:,1])) # confidence interval
#        # Elements for plot
#        xd = yd = [0,1]
#        q = self.ub - self.lb
#        #om = np.sum(c_table[:,0])/np.sum(c_table)   # overall (unconditional) relative frequency
#        #print(q,om)
#        #q =om ???? to be discussed
#        clim_x = clim_y = [q,q]
#        sk_line = [q/2,(1-q)/2+q]
#        #
#        fig = pl.figure(figsize=(7,5))
#        pl.plot(xd,yd,color='black',linestyle=':',linewidth=0.5)
#        pl.plot(xd,clim_y,color='black',linestyle=':',linewidth=0.5)
#        pl.plot(clim_x,yd,color='black',linestyle=':',linewidth=0.5)
#        pl.plot(xd,sk_line,color='black',linestyle='--',linewidth=0.5)
#        pl.fill_between(xd,xd,sk_line,facecolor='grey',alpha=0.2)
#        pl.fill_betweenx(yd,yd,clim_x,facecolor='grey',alpha=0.2)
#        pl.scatter(self.bins,p,s=np.sum(c_table,axis=1)/np.sum(c_table)*10000,color='deepskyblue',marker='o',alpha=0.5,edgecolors='none')
#        pl.plot(self.bins,p,color='deepskyblue',linestyle='-',linewidth=0.8,label='upper tercile')
#        pl.errorbar(self.bins,p,yerr=[p - ci_low,ci_upp - p],ecolor='deepskyblue',elinewidth=0.8,alpha=0.5)
#        pl.ylim(0.0,1.0)
#        pl.xlim(0.0,1.0)
#        pl.ylabel('Observed frequency \np(o|y)',fontsize=11)
#        pl.xlabel('Forecast probability \ny',fontsize=11)
#        #pl.legend(fontsize=7,loc='upper left')
#        pl.tight_layout()
#        pl.show()


        
        
