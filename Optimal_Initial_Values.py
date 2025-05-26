# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:17:52 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Relative_Errors import Relative_Errors
from Optimal_Tc import Optimal_Tc

def Optimal_Initial_Values (X, Y, M_tot, R_tot, L_tot, Tc):
    # We want to search in the interval R=[11.0,12.0] & L=[25.0,40.0]:
    dL = 10          
    dR = 0.5          

    L1 = L_tot-dL/2    # Minimun value for L_tot
    L2 = L_tot+dL      # Maximum value for L_tot

    R1 = R_tot-dR      # Minimun value for R_tot
    R2 = R_tot+dR      # Maximum value for R_tot
    
    nzooms = 5         # Number of zooms we want to make
    zoom = 1           # Number of the firts zoom         
    
    for n in range (0,nzooms):

        L_values = np.arange (L1, L2+(L2-L1)/10, (L2-L1)/10)  # Vector with the variance of L_tot:
        R_values = np.arange (R1, R2+(R2-R1)/10, (R2-R1)/10)  # Vector with the variance of R_tot:
        
        # Vectors with the total relative error at the boundary layer for each M_tot, L_tot:
        Erel_values = np.zeros ((len(L_values),len(R_values))) # Matrix for the relative errors
        Tc_values = np.zeros ((len(L_values),len(R_values)))  # Matrix for the optimal temperatures
        
        # Add values to each matrix element
        for n in range (0, len(L_values)):           # Rows 
            for m in range (0, len(R_values)):       # Columns   
                Tc_opt = Optimal_Tc (X, Y, M_tot, R_values[m], L_values[n], Tc)
                comparative_values = Relative_Errors (X, Y, M_tot, R_values[m], L_values[n], Tc_opt)
                Tc_values[n][m] = Tc_opt
                Erel_values[n][m] = comparative_values.at[2,'--']  
    
        # Look for the minimum error in the matrix:
        # Create a vector that includes the minimum values of the relative errors in each row:
        E_mins = []
        # Vector with the positions of the minimum errors in each row
        E_positions = []
        
        # Add the minimum row values and their positions:
        for n in range (0, len(L_values)):   
            E_mins.append(min(Erel_values[n,:]))
            E_positions.append(np.where(Erel_values[n,:] == E_mins[n])[0])
        
        # Set the minimum coordinates    
        L_min = int(np.where(E_mins == min(E_mins))[0])                  # Position of the L minimum
        R_min = int(E_positions[L_min])                                  # Position of the R minimum
        
        # Zoom Plot:
        ######################################################################
        plt.figure(figsize = (5, 5))

        # Show plot:
        plt.imshow(Erel_values, cmap = 'coolwarm')
        
        # Color bar
        cbar = plt.colorbar(label = "Relative Error (\%)")
        cbar.ax.tick_params(labelsize = 12)
        
        # Title & labels:
        plt.title(f'Zoom {zoom}', fontsize = 16, pad = 15)
        plt.xlabel(r'$R \ [10^{10} \ \mathrm{cm}]$', fontsize = 12)
        plt.ylabel(r'$L \ [10^{33} \ \mathrm{erg/s}]$', fontsize = 12)
        
        # Axis:
        plt.xticks(ticks = np.arange(len(R_values)), labels = [f"{val:.5f}" for val in R_values], rotation = -45, ha = 'left', fontsize = 11)
        plt.yticks(ticks = np.arange(len(L_values)), labels = [f"{val:.4f}" for val in L_values], fontsize = 11)
        
        # Tick axis:
        plt.tick_params(axis = 'both', which = 'major', length = 5, width = 1)
        plt.tick_params(bottom = True, top = True, left = True, right = True)
        
        # Invert Y axis:
        plt.gca().invert_yaxis()
        
        plt.show()
        ######################################################################
        
        # Set new values for new zoom:
        L_tot = L_values[L_min]
        dL = (L2-L1)/10
        R_tot = R_values[R_min]
        dR = (R2-R1)/10
        Tc = Tc_values[L_min][R_min]
        L1 = L_tot-dL
        L2 = L_tot+dL
        R1 = R_tot-dR
        R2 = R_tot+dR
        zoom = zoom+1
        
    return(R_tot, L_tot, Tc)