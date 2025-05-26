# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:16:55 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Relative_Errors import Relative_Errors

def Optimal_Tc (X, Y, M_tot, R_tot, L_tot, Tc):
    """
    The function calculates the optimal Tc for the star.
    
    Parameters:
    
    Returns:
    """
    zoom = Tc/10 # How big is the zoom

    # We are going to calculate 5 zooms:
    nzooms = 5
    for n in range (0, nzooms-1):
        loop = True
        while loop:
            Tc_values=np.arange(Tc-zoom, Tc+zoom, zoom/10) # Vector with the variance of Tc:

            #Vector with the total relative error at the boundary layer for each Tc:
            Erel_total=np.zeros(len(Tc_values))

            for i in range (0,len(Tc_values)):
                comparative_values = Relative_Errors (X, Y, M_tot, R_tot, L_tot, Tc_values[i])
                Erel_total[i] = comparative_values.at[2,'--']

            # Look for the minimum error:
            i_min = np.where (Erel_total == min(Erel_total))[0]

            # The minimum can't be at the extremes of the interval:
            if i_min == 0:
                zoom = zoom+zoom/10
            elif i_min == (len(Tc_values)-1):
                zoom = zoom+zoom/10
            else:
                loop = False

        Tc = Tc_values[i_min]
        Tc = Tc[0]
        zoom = zoom/10   
    return(Tc)