# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:18:29 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Energy_Production_Rate import Energy_Production_Rate


def Outer_Layers (X, Y, M_tot, R_tot, L_tot, Tc):
    E_exterior = []       # Energy production mechanism
    phase_exterior = []   # Model phase
    i_exterior = []       # Layer number 
    r_exterior = []       # Layer radius         [10^10 cm]
    P_exterior = []       # Layer pressure        [10^15 din cm^-2]
    T_exterior = []       # Layer temperature    [10^7 K]
    L_exterior = []       # Layer luminosity    [10^33 erg s^-1]
    M_exterior = []       # Layer mass           [10^33 g]
    
    fP_exterior = []    # f_i for pressure
    fT_exterior = []    # f_i for temperature
    fM_exterior = []    # f_i for mass
    fL_exterior = []    # f_i for luminosty
    
    nplus1_exterior = [] # n+1 value

    # Initial radius of the integration:
    R_ini = 0.9*R_tot   # (To avoid convergence problems) [10^10 cm]

    # Integration step:
    h = R_ini/100      # [10^10 cm]

    Z = 1-X-Y                    # Mass fraction of heavy elements                           [adimensional]
    mu = (1)/(2*X+3*Y/4+Z/2)     # Average molecular weight (constant throughout the star)   [adimensional]

    # Constants A1 & A2 necessary for the calculation of pressure and temperature in these first layers:
    A1 = 1.9022*mu*M_tot
    A2 = 10.645*np.sqrt(M_tot/(mu*Z*(X+1)*L_tot))
    
    # Consider mass & luminosity = constants
    L = L_tot
    M = M_tot 
    
    i = -1
    r = R_ini
    # Calculate the exterior layers:
    while r < R_tot:
        phase = '^^^^^^'
        r = R_ini+(-i)*h           
        T = A1*(1/r-1/R_tot)
        P = A2*(T**4.25)
        cycle, eps1, nu, X1, X2 = Energy_Production_Rate (T, P, X, Y)
    
        E_exterior.append(cycle)
        phase_exterior.append(phase)
        i_exterior.append(i)
        r_exterior.append(r)
        P_exterior.append(P)
        T_exterior.append(T)
        L_exterior.append(L)
        M_exterior.append(M)
        nplus1_exterior.append(0)

        # Constants Cm, Cp, Cl, Ct necesary to calculate the f_i's:
        Cm = 0.01523*mu
        Cp = 8.084*mu
        Cl = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)
        Ct = 0.01679*Z*(1+X)*(mu**2)
    
        # f_i factors:
        fP = -(Cp*P*M)/(T*(r**2))
        fT = -(Ct*(P**2)*L)/((T**8.5)*(r**2))
        fL = Cl*(P**2)*(T**(nu-2))*(r**2)
        fM = (Cm*P*(r**2))/(T)
    
        fP_exterior.append(fP)
        fT_exterior.append(fT)
        fL_exterior.append(fL)
        fM_exterior.append(fM)
        i = i-1
        
    
    exterior = pd.DataFrame ({
        'E': E_exterior,
        'Phase': phase_exterior,
        'i': i_exterior,
        'r': r_exterior,
        'P': P_exterior,
        'T': T_exterior,
        'L': L_exterior,
        'M': M_exterior,
    })
    
    return (exterior)