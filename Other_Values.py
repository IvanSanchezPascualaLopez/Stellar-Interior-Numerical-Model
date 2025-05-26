# -*- coding: utf-8 -*-
"""
Created on Thu May 22 19:44:45 2025

@author: sanch
"""

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

def Other_Values (X, Y, M_tot, R_tot, L_tot, Tc, star_values, boundary_values):
    # Effective Temperature Stimation [K]:
    T_eff = ((L_tot*(10**33)*(10**-7))/(4*np.pi*5.670374419e-8*(R_tot*10**8)**2))**(1/4)
    print('The effective temperature of the star, T_eff, is',int(T_eff),'K')
    ########################################################################
    # Emission peak (Wienn law) [nm]:
    lambda_max = (0.29/(T_eff))*1e7
    print('The emission peak of the stellar spectrum, lambda_max, is', round(lambda_max, 3), 'nm')
    #######################################################################
    # Superficial gravity:
    G = 6.6743e-8    # gravitational constant    [cm^3 g^-1 s^-2]  
    g = G*(M_tot*1e33)/((R_tot*1e10)**2)/100
    print ('The superficial gravity of the star, g, is', round(g, 3), 'm*s^-2')
    #######################################################################
    return T_eff, lambda_max, g