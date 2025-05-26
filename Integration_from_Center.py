# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:15:10 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Energy_Production_Rate import Energy_Production_Rate
from Integration_from_Surface import Integration_from_Surface


def Integration_from_Center (X, Y, M_tot, R_tot, L_tot, Tc):
    """
    The function calculates the convective nucleus of the star integrated from the center.
    This function also calculates the values at the radiative-convective boundary layer by applying a linear interpolation.
    Collects data for all layers in DataFrames.
    
    Parameters:
        X (float)
        Y (float)
        K (float)
        R_tot (float)
        Tc (float)
        boundary_from_surface (DataFrame)
    
    Returns:
        values_from_center (DataFrame): Includes next values:
            E (string):
            Phase (string):
            i (float):
            r (float):
            P (float):
            T (float):
            L (float):
            M (float):
            n+1 (float):
            
        boundary_from_center (DataFrame): Include next values:
            r boundary (float): Radius at the boundary layer [10^10 cm].
            P boundary (float): Pressure at the boundary layer [10^15 din cm^-2].
            T boundary (float): Temperature at the boundary layer [10^7 K].
            L boundary (float): Luminosity at the boundary layer [10^33 erg s^-1].
            M boundary (float): Mass at the boundary layer [10^33 g].
    """
    # Initial radius of the integration:
    R_ini = 0.9*R_tot   # (To avoid convergence problems) [10^10 cm]

    # Integration step:
    h = R_ini/100      # [10^10 cm]

    Z = 1-X-Y                    # Mass fraction of heavy elements                           [adimensional]
    mu = (1)/(2*X+3*Y/4+Z/2)     # Average molecular weight (constant throughout the star)   [adimensional]

    K, values_from_surface, boundary_from_surface = Integration_from_Surface (X, Y, M_tot, R_tot, L_tot)
    r_bound = boundary_from_surface.loc[0, 'r boundary']

    #####################################################################################################################################################
    
    E_prodc = []          # Energy Production Mechanism
    phasec = []           # Model phase
    i_valuesc = []       # layer number (begin in layer 0)
    r_valuesc = []       # layer radius          [10^10 cm]
    P_valuesc = []       # layer pressure        [10^15 din cm^-2]
    T_valuesc = []       # layer temperature    [10^7 K]
    L_valuesc = []       # layer luminosity    [10^33 erg s^-1]
    M_valuesc = []       # layer mass           [10^33 g]

    fT_valuesc = []    # f_i for temperature
    fM_valuesc = []    # f_i for mass
    fL_valuesc = []    # f_i for luminosity

    nplus1_valuesc = []

    # First 3  layers:
    layers = 3 # Number of layers
    
    # Calculate the values in each layer:
    for i in range (0, layers):
        phase = 'CENTER'
        r = i*h
        M = 0.005077*mu*K*(Tc**1.5)*(r**3)
        T = Tc-0.008207*(mu**2)*K*(Tc**1.5)*(r**2)
        P = K*(T**2.5)
        cycle, eps1, nu, X1, X2 = Energy_Production_Rate (T, P, X, Y)   # Energy production parameters
        L = 0.006150*eps1*X1*X2*(10**nu)*(mu**2)*(K**2)*(Tc**(3+nu))*(r**3)

        E_prodc.append('--')
        phasec.append(phase)
        i_valuesc.append(i)
        r_valuesc.append(r)
        P_valuesc.append(P)
        T_valuesc.append(T)
        L_valuesc.append(L)
        M_valuesc.append(M)
        nplus1_valuesc.append(0)

        # Constants Cm, Cp, Cl, Ct necessaries for the calculos of the f_i's:
        Cm = 0.01523*mu
        Cl = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)
        Ct2 = 3.234*mu

        # f_i factors:
        fT = -(Ct2*M)/(r**2)
        fL = Cl*(K**2)*(T**(3+nu))*(r**2)
        fM = Cm*K*(T**1.5)*(r**2)

        fT_valuesc.append(fT)
        fL_valuesc.append(fL)
        fM_valuesc.append(fM)
    
    #####################################################################################################################################################

    # Integration from the center for the convective nucleus:

    # We have already calculated the 3 initial layers (i = 0, 1 & 2)
    i = 2 # Number of the last calculated layer

    #We use the logical variables loop1, loop2 & loop3 to control each of the loops needed to program the algorithm. (convective nucleus):

    loop1 = True
    while loop1:
        # Executing Step 1
        r = r_valuesc[i]+h           # radius of the layer i+1
        # Executing Step 2bis
        Deltai1T = h*(fT_valuesc[i]-fT_valuesc[i-1])
        Test = T_valuesc[i]+h*fT_valuesc[i]+0.5*Deltai1T 
    
        loop2 = True
        while loop2:       
            # Executing Politrope Step
            Pest = K*(Test**2.5)

            # Executing Step 3
            fM = 0.01523*mu*(Pest/Test)*(r**2)
            Deltaii1M = h*(fM-fM_valuesc[i])
            Mcal = M_valuesc[i]+h*fM-0.5*Deltaii1M
        
            # Executing Step 7bis
            if r > 0:
                fT = -3.234*mu*(Mcal/(r**2))
                Deltaii1T = h*(fT-fT_valuesc[i])
                Tcal = T_valuesc[i]+h*fT-0.5*Deltaii1T
            else:
                Tcal = Test
        
            # Executing Step 8
            ErelT = abs(Tcal-Test)/Tcal
            Erelmax = 0.0001
            if ErelT < Erelmax:
                loop2 = False
            else:
                Test = Tcal
            
        # executing Politrope Step
        Pcal = K*(Tcal**2.5)
    
        # Executing Step 6bis
        cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Tcal, Pcal, X, Y)
        fL = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)*(Pcal**2)*(Tcal**(nu-2))*(r**2)
        Deltaii1L = h*(fL-fL_valuesc[i])
        Deltaii2L = Deltaii1L-h*(fL_valuesc[i]-fL_valuesc[i-1])
        Lcal = L_valuesc[i]+h*fL-0.5*Deltaii1L-(1/12)*Deltaii2L

    
        if r >= r_bound:
            boundary = i+1
            loop1 = False
        else:
            phase = 'CONV'
            cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Tcal, Pcal, X, Y)

            E_prodc.append(cycle)
            phasec.append(phase)
            i_valuesc.append(i+1)
            r_valuesc.append(r)
            P_valuesc.append(Pcal)
            T_valuesc.append(Tcal)
            L_valuesc.append(Lcal)
            M_valuesc.append(Mcal)
            nplus1_valuesc.append(0)

            fT_valuesc.append(fT)
            fL_valuesc.append(fL)
            fM_valuesc.append(fM)
        
            i += 1   # Next layer
            #Iterate loop1 withnext layer"


    phase = 'CONV'
    cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Tcal, Pcal, X, Y)

    E_prodc.append(cycle)
    phasec.append(phase)
    i_valuesc.append(i+1)
    r_valuesc.append(r)
    P_valuesc.append(Pcal)
    T_valuesc.append(Tcal)
    L_valuesc.append(Lcal)
    M_valuesc.append(Mcal)
    nplus1_valuesc.append(0)

    fT_valuesc.append(fT)
    fL_valuesc.append(fL)
    fM_valuesc.append(fM)

    #####################################################################################################################################################

    # Boundary values (linear interpolation):
    r_boundaryc = []
    P_boundaryc = []
    T_boundaryc = []
    L_boundaryc = []
    M_boundaryc = []

    r_boundc = r_bound
    P_boundc = ((P_valuesc[boundary]-P_valuesc[boundary-1])/(r_valuesc[boundary]-r_valuesc[boundary-1]))*(r_boundc-r_valuesc[boundary-1])+P_valuesc[boundary-1]
    T_boundc = ((T_valuesc[boundary]-T_valuesc[boundary-1])/(r_valuesc[boundary]-r_valuesc[boundary-1]))*(r_boundc-r_valuesc[boundary-1])+T_valuesc[boundary-1]
    L_boundc = ((L_valuesc[boundary]-L_valuesc[boundary-1])/(r_valuesc[boundary]-r_valuesc[boundary-1]))*(r_boundc-r_valuesc[boundary-1])+L_valuesc[boundary-1]
    M_boundc = ((M_valuesc[boundary]-M_valuesc[boundary-1])/(r_valuesc[boundary]-r_valuesc[boundary-1]))*(r_boundc-r_valuesc[boundary-1])+M_valuesc[boundary-1]

    r_boundaryc.append(r_boundc)
    P_boundaryc.append(P_boundc)
    T_boundaryc.append(T_boundc)
    L_boundaryc.append(L_boundc)
    M_boundaryc.append(M_boundc)
    
    #####################################################################################################################################################

    values_from_center = pd.DataFrame({
        'E': E_prodc,
        'Phase': phasec,
        'i': i_valuesc,
        'r': r_valuesc,
        'P': P_valuesc,
        'T': T_valuesc,
        'L': L_valuesc,
        'M': M_valuesc,
        'n+1': nplus1_valuesc,
    })

    boundary_from_center = pd.DataFrame ({
        'i boundary': boundary-1,
        'r boundary': r_boundaryc,
        'P boundary': P_boundaryc,
        'T boundary': T_boundaryc,
        'L boundary': L_boundaryc,
        'M boundary': M_boundaryc,   
    })

    #####################################################################################################################################################

    return (values_from_surface, boundary_from_surface, values_from_center, boundary_from_center)