# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:16:01 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Integration_from_Center import Integration_from_Center

def Relative_Errors (X, Y, M_tot, R_tot, L_tot, Tc):
    """
    The function calculates the optimal Tc for the star.
    
    Parameters:
    
    Returns:
    """

    values_from_surface, boundary_from_surface, values_from_center, boundary_from_center = Integration_from_Center (X, Y, M_tot, R_tot, L_tot, Tc)
    
    column1_comparative = []
    column2_comparative = []
    r_comparative = []
    P_comparative = []
    T_comparative = []
    L_comparative = []
    M_comparative = []
    column8_comparative = []
    column9_comparative = []

    column1_comparative.append('DOWN')
    column2_comparative.append('->')
    r_comparative.append(boundary_from_surface['r boundary'].iloc[0])
    P_comparative.append(boundary_from_surface['P boundary'].iloc[0])
    T_comparative.append(boundary_from_surface['T boundary'].iloc[0])
    L_comparative.append(boundary_from_surface['L boundary'].iloc[0])
    M_comparative.append(boundary_from_surface['M boundary'].iloc[0])
    column8_comparative.append('')
    column9_comparative.append('')

    column1_comparative.append('UP')
    column2_comparative.append('->')
    r_comparative.append(boundary_from_center['r boundary'].iloc[0])
    P_comparative.append(boundary_from_center['P boundary'].iloc[0])
    T_comparative.append(boundary_from_center['T boundary'].iloc[0])
    L_comparative.append(boundary_from_center['L boundary'].iloc[0])
    M_comparative.append(boundary_from_center['M boundary'].iloc[0])
    column8_comparative.append('')
    column9_comparative.append('')

    ErelP = (abs(boundary_from_surface['P boundary'].iloc[0]-boundary_from_center['P boundary'].iloc[0])/boundary_from_surface['P boundary'].iloc[0])*100
    ErelT = (abs(boundary_from_surface['T boundary'].iloc[0]-boundary_from_center['T boundary'].iloc[0])/boundary_from_surface['T boundary'].iloc[0])*100
    ErelL = (abs(boundary_from_surface['L boundary'].iloc[0]-boundary_from_center['L boundary'].iloc[0])/boundary_from_surface['L boundary'].iloc[0])*100
    ErelM = (abs(boundary_from_surface['M boundary'].iloc[0]-boundary_from_center['M boundary'].iloc[0])/boundary_from_surface['M boundary'].iloc[0])*100
    Erel_totp = np.sqrt((ErelP**2)+(ErelT**2)+(ErelL**2)+(ErelM**2))

    column1_comparative.append('rel.error')
    column2_comparative.append('(%):')
    r_comparative.append('')
    P_comparative.append(ErelP)
    T_comparative.append(ErelT)
    L_comparative.append(ErelL)
    M_comparative.append(ErelM)
    column8_comparative.append('TOTAL.REL.ERROR:')
    column9_comparative.append(Erel_totp)

    comparative_values = pd.DataFrame ({
        '': column1_comparative,
        '-': column2_comparative,
        'r': r_comparative,
        'P': P_comparative,
        'T': T_comparative,
        'L': L_comparative,
        'M': M_comparative, 
        '---': column8_comparative,
        '--': column9_comparative
    })

    return(comparative_values)
