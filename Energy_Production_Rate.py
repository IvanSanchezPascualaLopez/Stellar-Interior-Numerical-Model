# -*- coding: utf-8 -*-
"""
Created on Tue May  6 19:39:18 2025

@author: sanch
"""


# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

def Energy_Production_Rate (T,P,X,Y):     
    """
    The function calculates the rate of energy production inside of the star based on local conditions.
    It decides the predominant energy production cycle (pp or CNO).
    
    Parameters:
        T (float): Temperature [10^7 K].
        P (float): Pressure [10^15 din cm^-2].
        X (float): Mass fraction of hydrogen [adimensional].
        Y (float): Mass fraction of helium [adimensional].

    Returns:
        cycle (string): shows the predominant cycle, "pp" or "CNO". Also shows if there is no energy production, "--".
        eps1, nu, X1, X2 (float): energy production rate parameters.
    """
    Z = 1-X-Y                    # Mass fraction of heavy elements                           [adimensional]
    mu = (1)/(2*X+3*Y/4+Z/2)     # Average molecular weight (constant throughout the star)   [adimensional]
    Na = 6.022e23                # Avogadro number                                           [adimensional]
    H = 1/Na                     # Avogadro numer inverse                                    [adimensional]
    k = 1.380649e-23             # Boltzmann constant                                        [din*cm/10^7K]=[J/K]
    rho = ((mu*H*P)/(k*T))       # Density                                                   [g/cm^3]
    
    # Parameters of the energy production cycles:
    X1_pp=X
    X2_pp=X
    X1_CNO=X
    X2_CNO=Z/3

    # nu and epsilon_1 values as a function of temperature for the pp cycle:
    if T < 0.4:
        eps1_pp = 0
        nu_pp = 0
    elif T >= 0.4 and T < 0.6:
        eps1_pp = 10**(-6.84)
        nu_pp = 6
    elif T >= 0.6 and T<0.95:
        eps1_pp = 10**(-6.04)
        nu_pp = 5
    elif T >= 0.95 and T < 1.2:
        eps1_pp = 10**(-5.56)
        nu_pp = 4.5
    elif T >= 1.2 and T<1.65:
        eps1_pp = 10**(-5.02)
        nu_pp = 4
    elif T >= 1.65 and T < 2.4:
        eps1_pp = 10**(-4.40)
        nu_pp = 3.5
    else:
        eps1_pp = 0
        nu_pp = 0
        
    # nu and epsilon_1 values as a function of temperature for the CNO cycle:
    if T < 1.2:
        eps1_CNO = 0
        nu_CNO = 0
    elif T >= 1.2 and T < 1.6:
        eps1_CNO = 10**(-22.2)
        nu_CNO = 20
    elif T >= 1.6 and T < 2.25:
        eps1_CNO = 10**(-19.8)
        nu_CNO = 18
    elif T >= 2.25 and T<2.75:
        eps1_CNO = 10**(-17.1)
        nu_CNO = 16
    elif T >= 2.75 and T < 3.6:
        eps1_CNO = 10**(-15.6)
        nu_CNO = 15
    elif T >= 3.6 and T < 5:
        eps1_CNO = 10**(-12.5)
        nu_CNO = 13
    else:
        eps1_CNO = 0
        nu_CNO = 0

    # Energy cycles formulas:
    pp = eps1_pp*X1_pp*X2_pp*rho*(10*T)**nu_pp
    CNO = eps1_CNO*X1_CNO*X2_CNO*rho*(10*T)**nu_CNO

    # Compare both values and adopt the one corresponding to the cycle that produces the most energy:
    if pp > CNO:
        cycle = 'pp'
        eps1 = eps1_pp
        nu = nu_pp
        X1 = X1_pp
        X2 = X2_pp
    else:
        cycle = 'CN'
        eps1 = eps1_CNO
        nu = nu_CNO
        X1 = X1_CNO
        X2 = X2_CNO
    if eps1 == 0:
        cycle = '--'
        
    return cycle, eps1, nu, X1, X2 
