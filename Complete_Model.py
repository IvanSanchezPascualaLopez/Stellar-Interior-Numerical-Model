# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:19:10 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Integration_from_Center import Integration_from_Center
from Outer_Layers import Outer_Layers
pd.set_option ('display.max_rows', None)

def Complete_Model (X, Y, M_tot, R_tot, L_tot, Tc):
    R_tot
    E_complete = []    # Energy Production Mechanism
    phase_complete = [] # Model phase
    i_complete = []       # Layer number 
    r_complete = []       # Layer radius          [10^10 cm]
    P_complete = []       # layer pressure       [10^15 din cm^-2]
    T_complete = []       # layer temperature    [10^7 K]
    L_complete = []       # layer luminosity   [10^33 erg s^-1]
    M_complete = []       # layer mass           [10^33 g]
    
    fP_complete = []    # f_i for pressure
    fT_complete = []    # f_i for temperature
    fM_complete = []    # f_i for mass
    fL_complete = []    # f_i for luminosity
    
    nplus1_complete = []

    # Add the outer layers:
    outer = Outer_Layers (X, Y, M_tot, R_tot, L_tot, Tc)
    n = len(outer['E'])-2
    while n >= 0:
        E_complete.append(outer.at[n,'E'])
        phase_complete.append(outer.at[n,'Phase'])
        i_complete.append(outer.at[n,'i'])
        r_complete.append(outer.at[n,'r'])
        P_complete.append(outer.at[n,'P'])
        T_complete.append(outer.at[n,'T'])
        L_complete.append(outer.at[n,'L'])
        M_complete.append(outer.at[n,'M'])
        nplus1_complete.append(0)
    
        n = n-1

    # Add the radiative layers:
    values_from_surface, boundary_from_surface, values_from_center, boundary_from_center = Integration_from_Center (X, Y, M_tot, R_tot, L_tot, Tc)
    n = 0
    while n <= boundary_from_surface.at[0,'i boundary']:
        E_complete.append(values_from_surface.at[n,'E'])
        phase_complete.append(values_from_surface.at[n,'Phase'])
        i_complete.append(values_from_surface.at[n,'i'])
        r_complete.append(values_from_surface.at[n,'r'])
        P_complete.append(values_from_surface.at[n,'P'])
        T_complete.append(values_from_surface.at[n,'T'])
        L_complete.append(values_from_surface.at[n,'L'])
        M_complete.append(values_from_surface.at[n,'M'])
        nplus1_complete.append(values_from_surface.at[n,'n+1'])
    
        n = n+1

    # Add the convective layers:
    i = boundary_from_surface.at[0,'i boundary']+1
    n = boundary_from_center.at[0,'i boundary']
    while n >= 0:
        E_complete.append(values_from_center.at[n,'E'])
        phase_complete.append(values_from_center.at[n,'Phase'])
        i_complete.append(i)
        r_complete.append(values_from_center.at[n,'r'])
        P_complete.append(values_from_center.at[n,'P'])
        T_complete.append(values_from_center.at[n,'T'])
        L_complete.append(values_from_center.at[n,'L'])
        M_complete.append(values_from_center.at[n,'M'])
        nplus1_complete.append(values_from_center.at[n,'n+1'])
    
        i = i+1
        n = n-1
    
    star_values = pd.DataFrame ({
        'E': E_complete,
        'Fase': phase_complete,
        'i': i_complete,
        'r': r_complete,
        'P': P_complete,
        'T': T_complete,
        'L': L_complete,
        'M': M_complete,
        'n+1': nplus1_complete
    })
    print (star_values)
    return star_values, boundary_from_surface