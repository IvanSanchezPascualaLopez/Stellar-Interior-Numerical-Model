# -*- coding: utf-8 -*-
"""
Created on Tue May  6 19:36:51 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Initial_Values import Initial_Values
from Energy_Production_Rate import Energy_Production_Rate
from Integration_from_Surface import Integration_from_Surface
from Integration_from_Center import Integration_from_Center
from Relative_Errors import Relative_Errors
from Optimal_Initial_Values import Optimal_Initial_Values
from Optimal_Tc import Optimal_Tc
from Outer_Layers import Outer_Layers
from Complete_Model import Complete_Model
from Parameters_vs_Radius import Parameters_vs_Radius
from Other_Values import Other_Values
from Hipparcos_HR_Diagram import Hipparcos_HR_Diagram

# Use LaTeX font:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",      
    "font.serif": ["Computer Modern Roman"],  
})

pd.set_option('display.max_rows',None)      # Display full DataFrame rows
pd.set_option('display.max_columns',None)   # Display full DataFrame columns
###############################################################################
# Initial parameters (Test Model):
#M_tot=5.0       # Total mass of the star        [10^33 g]    
#X=0.75          # Mass fraction of hydrogen     [adimensional]
#Y=0.22          # Mass fraction of helium       [adimensional]

# Initial values (Test Model):
#R_tot=11.5      # Total radius                  [10^10 cm] 
#L_tot=70        # Total luminosity              [10^33 erg s^-1]
#Tc=2.0          # Central temperature           [10^7 K]
############################################################################
# Initial parameters (My Star):
#M_tot=4.8       # Total mass of the star        [10^33 g]    
#X=0.75          # Mass fraction of hydrogen     [adimensional]
#Y=0.20          # Mass fraction of helium       [adimensional]

# Initial values (My Star):
#R_tot=11.5      # Total radius                  [10^10 cm] 
#L_tot=30.0        # Total luminosity            [10^33 erg s^-1]
#Tc=1.5          # Central temperature           [10^7 K]
############################################################################
# Choose the initial values:
X, Y, M_tot, R_tot, L_tot, Tc = Initial_Values ()  
##########################################################################
# Look for the optimal initial values:
R_tot, L_tot, Tc = Optimal_Initial_Values (X, Y, M_tot, R_tot, L_tot, Tc) 
###############################################################################
print('The optimal value of L_tot is :',L_tot,'·10^33 erg s^-1')
print('The optimal value of R_tot is :',R_tot,'·10^10 cm')
print('The optimal value of Tc is :',Tc,'·10^7 K')
###########################################################################
# Calculate the relative errors:
comparative_values = Relative_Errors (X, Y, M_tot, R_tot, L_tot, Tc) 
print(comparative_values)
###########################################################################
# Calculate the layer values of the star
star_values, boundary_values = Complete_Model (X, Y, M_tot, R_tot, L_tot, Tc) 
print(boundary_values)
#############################################################################
# Show the parameter dependencies with radius:
rho, eps_values, kappa = Parameters_vs_Radius (X, Y, M_tot, R_tot, L_tot, Tc, star_values, boundary_values)
##########################################################################
# Calculate other relevant values of the star:
T_eff, lambda_max, g = Other_Values (X, Y, M_tot, R_tot, L_tot, Tc, star_values, boundary_values)
###########################################################################
# Show the star on a HR Diagram (Hipparcos)
Hipparcos_HR_Diagram (T_eff, L_tot)
