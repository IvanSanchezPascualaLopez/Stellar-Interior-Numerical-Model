# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:14:31 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Energy_Production_Rate import Energy_Production_Rate

def Integration_from_Surface (X, Y, M_tot, R_tot, L_tot):
    """
    The function calculates the radiative envelope and convective nucleus of the star integrated from the surface.
    This function also calculates the values at the radiative-convective boundary layer by applying a linear interpolation.
    Collects data for all layers in DataFrames.
    
    Parameters:
        X (float): Mass fraction of hydrogen [adimensional].
        Y (float): Mass fraction of helium [adimensional].
        M_tot (float): Total mass of the star [10^33 g]
        R_tot (float): Total radius of the star [10^10 cm].
        L_tot (float): Total luminosity of the star [10^33 erg s^-1].

    Returns:
        K (float): Polytropic Constant. Necessary for the integration from the center.
        values_from_surface (DataFrame): Includes next values:
            E (string):
            Phase (string):
            i (float):
            r (float):
            P (float):
            T (float):
            L (float):
            M (float):
            n+1 (float):
            
        boundary_from_surface (DataFrame): Include next values:
            r boundary (float): Radius at the boundary layer [10^10 cm].
            P boundary (float): Pressure at the boundary layer [10^15 din cm^-2].
            T boundary (float): Temperature at the boundary layer [10^7 K].
            L boundary (float): Luminosity at the boundary layer [10^33 erg s^-1].
            M boundary (float): Mass at the boundary layer [10^33 g].
    """
    # Initial radius of the integration:
    R_ini = 0.9*R_tot   # (To avoid convergence problems) [10^10 cm]

    # Integration step:
    h =- R_ini/100      # [10^10 cm]

    Z = 1-X-Y                    # Mass fraction of heavy elements                           [adimensional]
    mu = (1)/(2*X+3*Y/4+Z/2)     # Average molecular weight (constant throughout the star)   [adimensional]

    #####################################################################################################################################################
    
    # Create the lists to store the values of each layer:
    E_prod = []         # Energy Production Mechanism
    model_phase = []    # Model phase
    i_values = []       # layer numbers 
    r_values = []       # Layer radius           [10^10 cm]
    P_values = []       # Layer pressures        [10^15 din cm^-2]
    T_values = []       # Layer temperatures     [10^7 K]
    L_values = []       # layer luminosities     [10^33 erg s^-1]
    M_values = []       # Layer mass             [10^33 g]
    nplus1_values = []  # n+1 values             [adimensional]

    # Create the lists to store the f_i values of each layer:
    fP_values=[]    # f_i for pressure
    fT_values=[]    # f_i for temperature
    fM_values=[]    # f_i for mass
    fL_values=[]    # f_i for luminosity

    #####################################################################################################################################################

    # First 3 superficial layers (M,l=cons; fM=fL=0):
    layers = 3 # Number of layers

    # Constants A1 & A2 necessary for the calculation of pressure and temperature in these first layers:
    A1 = 1.9022*mu*M_tot
    A2 = 10.645*np.sqrt(M_tot/(mu*Z*(X+1)*L_tot))

    # Luminosity and mass are constant and equal to the total mass and luminosity in these first three layers:
    L = L_tot
    M = M_tot 

    # Calculate the values in each layer:
    for i in range (0, layers):
        phase = 'BEGIN'
        r = R_ini+(i)*h
        T = A1*(1/r-1/R_tot)
        P = A2*(T**4.25)
        cycle, eps1 , nu, X1, X2 = Energy_Production_Rate (T, P, X, Y)

        E_prod.append(cycle)
        model_phase.append(phase)
        i_values.append(i)
        r_values.append(r)
        P_values.append(P)
        T_values.append(T)
        L_values.append(L)
        M_values.append(M)
        nplus1_values.append(0)

        # Constants Cm, Cp, Cl, Ct necesary to calculate the f_i's:
        Cm = 0                              # Mass does not vary
        Cp = 8.084*mu
        Cl = 0                              # Luminosity does not vary
        Ct = 0.01679*Z*(1+X)*(mu**2)

        # f_i factors:
        fP = -(Cp*P*M)/(T*(r**2))
        fT = -(Ct*(P**2)*L)/((T**8.5)*(r**2))
        fL = Cl*(P**2)*(T**(nu-2))*(r**2)
        fM = (Cm*P*(r**2))/(T)

        fP_values.append(fP)
        fT_values.append(fT)
        fL_values.append(fL)
        fM_values.append(fM)

    #####################################################################################################################################################

    # Integration from the surface for the radiative envelope:

    # We have already calculated the 3 initial layers (i = 0, 1 & 2)
    i = 2 # Number of the last calculated layer

    #We use the logical variables loop1, loop2 & loop3 to control each of the loops needed to program the algorithm of phase A.1. (radiative envelope):
    loop1 = True
    while loop1:
        # Executing Step 1
        r = r_values[i]+h           # Radius for the i+1 layer
    
        # Executing Step 2
        # Calculate the differences:
        Deltai1P = h*(fP_values[i]-fP_values[i-1])
        Deltai2P = Deltai1P-h*(fP_values[i-1]-fP_values[i-2])
        Deltai1T = h*(fT_values[i]-fT_values[i-1])
        
        # Calculate the estimated pressure and temperature for the i+1 layer:
        Pest = P_values[i]+h*fP_values[i]+0.5*Deltai1P+(5/12)*Deltai2P
        Test = T_values[i]+h*fT_values[i]+0.5*Deltai1T
    
        loop2 = True
        while loop2:
            loop3 = True
            while loop3:
                # Executing Step 3
                fM = 0.01523*mu*(Pest/Test)*(r**2)
                Deltaii1M = h*(fM-fM_values[i])
                Mcal = M_values[i]+h*fM-0.5*Deltaii1M   # Calculated mass 
            
                # Executing Step 4
                fP = -8.084*mu*(Pest/Test)*(Mcal/(r**2))
                Deltaii1P = h*(fP-fP_values[i])
                Pcal = P_values[i]+h*fP-0.5*Deltaii1P    # Calculated pressure
            
                # Executing Step 5
                ErelP = abs(Pcal-Pest)/Pcal
                Erelmax = 0.0001
                if ErelP < Erelmax:
                    loop3 = False
                else:
                    Pest = Pcal
        
            # Executing Step 6
            cycle, eps1 , nu, X1, X2 = Energy_Production_Rate (Test, Pcal, X, Y)   # Energy production parameters
            fL = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)*(Pcal**2)*(Test**(nu-2))*(r**2)
            Deltaii1L = h*(fL-fL_values[i])
            Deltaii2L = Deltaii1L-h*(fL_values[i]-fL_values[i-1])
            Lcal = L_values[i]+h*fL-0.5*Deltaii1L-(1/12)*Deltaii2L  # Calculated luminosity
        
            # Executing Step 7
            fT = -0.01679*Z*(1+X)*(mu**2)*((Pcal**2)/(Test**8.5))*(Lcal/(r**2))
            Deltaii1T = h*(fT-fT_values[i])
            Tcal = T_values[i]+h*fT-0.5*Deltaii1T   # Calcculated temperature
        
            # Executing Step 8
            ErelT = abs(Tcal-Test)/Tcal
            if ErelT < Erelmax:
                loop2 = False
            else:
                Test = Tcal
                Pest = Pcal
            
        # Executing Step 9
        nplus1 = (Tcal/Pcal)*((h*fP)/(h*fT))
    
        # Executing Step 10
        if nplus1 <= 2.5:
            boundary = i # Radiative-convective boundary layer (this layer is not radiative).
            loop1 = False
        else:
            phase = 'RADIA'
            cycle,eps1,nu, X1, X2 = Energy_Production_Rate (Test,Pcal,X,Y)

            E_prod.append(cycle)
            model_phase.append(phase)
            i_values.append(i+1)
            r_values.append(r)
            P_values.append(Pcal)
            T_values.append(Tcal)
            L_values.append(Lcal)
            M_values.append(Mcal)
            nplus1_values.append(nplus1)

            fP_values.append(fP)
            fT_values.append(fT)
            fL_values.append(fL)
            fM_values.append(fM)
        
            i += 1   # Next layer
            #Iterate loop1 for the next layer"

    #####################################################################################################################################################
    
    
    # Executing Step 2
    # Calculate the differences:
    Deltai1P = h*(fP_values[boundary]-fP_values[boundary-1])
    Deltai2P = Deltai1P-h*(fP_values[boundary-1]-fP_values[boundary-2])
    Deltai1T = h*(fT_values[boundary]-fT_values[boundary-1])
    
    # Calculate the pressure and temperature for the convenctive layer:
    Pest = P_values[boundary]+h*fP_values[boundary]+0.5*Deltai1P+(5/12)*Deltai2P
    Test = T_values[boundary]+h*fT_values[boundary]+0.5*Deltai1T

    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
            # Executing Step 3
            fM = 0.01523*mu*(Pest/Test)*(r**2)
            Deltaii1M = h*(fM-fM_values[boundary])
            Mcal = M_values[boundary]+h*fM-0.5*Deltaii1M
            
            # Executing Step 4
            fP = -8.084*mu*(Pest/Test)*(Mcal/(r**2))
            Deltaii1P = h*(fP-fP_values[boundary])
            Pcal = P_values[boundary]+h*fP-0.5*Deltaii1P
            
            # Executing Step 5
            ErelP = abs(Pcal-Pest)/Pcal
            Erelmax = 0.0001
            if ErelP < Erelmax:
                loop3 = False
            else:
                Pest=Pcal
        
        # Executing Step 6
        cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Test, Pcal, X, Y)   # Energy production parameters
        fL = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)*(Pcal**2)*(Test**(nu-2))*(r**2)
        Deltaii1L = h*(fL-fL_values[boundary])
        Deltaii2L = Deltaii1L-h*(fL_values[boundary]-fL_values[boundary-1])
        Lcal = L_values[boundary]+h*fL-0.5*Deltaii1L-(1/12)*Deltaii2L
        
        # Executing Step 7
        fT = -0.01679*Z*(1+X)*(mu**2)*((Pcal**2)/(Test**8.5))*(Lcal/(r**2))
        Deltaii1T = h*(fT-fT_values[boundary])
        Tcal = T_values[boundary]+h*fT-0.5*Deltaii1T
        
        # Executing Step 8
        ErelT = abs(Tcal-Test)/Tcal
        if ErelT < Erelmax:
            loop2 = False
        else:
            Test = Tcal
            Pest = Pcal
        
    K = Pcal/(Tcal**2.5)  # Will be constant for every layer

    #####################################################################################################################################################

    # Already calculated the radiative envelope.

    i = boundary  # Numer of the last calculated radiative layer. 

    # We use the logical variables loop1, loop2 & loop3 to control each of the loops needed to program the algorithm of phase A.2. (convective nucleus):

    loop1 = True
    while loop1:
        # Executing Step 1
        r = r_values[i]+h           # Radius for the i+1 layer
        # Executing Step 2bis
        Deltai1T = h*(fT_values[i]-fT_values[i-1])
        Test = T_values[i]+h*fT_values[i]+0.5*Deltai1T 
    
        loop2 = True
        while loop2:       
            # Executing Politrope Step
            if Test <= 0.000001 or np.isnan(Test) or np.isinf(Test): 
                break
            Pest = K*(Test**2.5)

            # Executing Step 3
            fM = 0.01523*mu*(Pest/Test)*(r**2)
            Deltaii1M = h*(fM-fM_values[i])
            Mcal = M_values[i]+h*fM-0.5*Deltaii1M
        
            # Executing Step 7bis
            if r > 0.000001:
                fT = -3.234*mu*(Mcal/(r**2))
                Deltaii1T = h*(fT-fT_values[i])
                Tcal = T_values[i]+h*fT-0.5*Deltaii1T
                if Tcal <= 0.000001 or np.isnan(Tcal) or np.isinf(Tcal):
                    break
            else:
                Tcal=Test
        
            # Executing Step 8
            ErelT = abs(Tcal-Test)/Tcal
            if ErelT < Erelmax:
                loop2 = False
            else:
                Test = Tcal
            
        # Executing Politrope Step
        Pcal = K*(Tcal**2.5)
    
        #Executing Step 6bis
        cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Tcal, Pcal, X, Y)   # Energy production parameters
        fL = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)*(Pcal**2)*(Tcal**(nu-2))*(r**2)
        Deltaii1L = h*(fL-fL_values[i])
        Deltaii2L = Deltaii1L-h*(fL_values[i]-fL_values[i-1])
        Lcal = L_values[i]+h*fL-0.5*Deltaii1L-(1/12)*Deltaii2L

    
        if r < 0.000001:
            loop1 = False
        else:
            phase = 'CONV'
            cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Tcal, Pcal, X, Y)

            E_prod.append(cycle)
            model_phase.append(phase)
            i_values.append(i+1)
            r_values.append(r)
            P_values.append(Pcal)
            T_values.append(Tcal)
            L_values.append(Lcal)
            M_values.append(Mcal)
            nplus1_values.append(0)

            fP_values.append(fP)
            fT_values.append(fT)
            fL_values.append(fL)
            fM_values.append(fM)
        
            i += 1   # Next layer
            #Iterate loop1 for the new layer"

    phase = 'CONV'
    cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Tcal, Pcal, X, Y)

    E_prod.append(cycle)
    model_phase.append(phase)
    i_values.append(i+1)
    r_values.append(0)
    P_values.append(Pcal)
    T_values.append(Tcal)
    L_values.append(Lcal)
    M_values.append(Mcal)
    nplus1_values.append(0)

    fP_values.append(fP)
    fT_values.append(fT)
    fL_values.append(fL)
    fM_values.append(fM)

    #####################################################################################################################################################

    # n+1 value for the first convective layer:
    # Excecuting Step 1
    r = r_values[boundary]+h           # Radius of layer i+1
    
    # Executing Step 2
    # Calculate the differences:
    Deltai1P = h*(fP_values[boundary]-fP_values[boundary-1])
    Deltai2P = Deltai1P-h*(fP_values[boundary-1]-fP_values[boundary-2])
    Deltai1T = h*(fT_values[boundary]-fT_values[boundary-1])

    # Calculate pressure and temperature for layer i+1:
    Pest = P_values[boundary]+h*fP_values[boundary]+0.5*Deltai1P+(5/12)*Deltai2P
    Test = T_values[boundary]+h*fT_values[boundary]+0.5*Deltai1T
    
    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
            # Executing Step 3
            fM = 0.01523*mu*(Pest/Test)*(r**2)
            Deltaii1M = h*(fM-fM_values[boundary])
            Mcal = M_values[boundary]+h*fM-0.5*Deltaii1M
            
            # Executing Step 4
            fP = -8.084*mu*(Pest/Test)*(Mcal/(r**2))
            Deltaii1P = h*(fP-fP_values[boundary])
            Pcal = P_values[boundary]+h*fP-0.5*Deltaii1P
            
            # Executing Step 5
            ErelP = abs(Pcal-Pest)/Pcal
            Erelmax = 0.0001
            if ErelP < Erelmax:
                loop3 = False
            else:
                Pest = Pcal
        
        # Executing Step 6
        cycle, eps1, nu, X1, X2 = Energy_Production_Rate (Test, Pcal, X, Y)
        fL = 0.01845*eps1*X1*X2*(10**nu)*(mu**2)*(Pcal**2)*(Test**(nu-2))*(r**2)
        Deltaii1L = h*(fL-fL_values[boundary])
        Deltaii2L = Deltaii1L-h*(fL_values[boundary]-fL_values[boundary-1])
        Lcal = L_values[boundary]+h*fL-0.5*Deltaii1L-(1/12)*Deltaii2L
        
        # Executing Step 7
        fT = -0.01679*Z*(1+X)*(mu**2)*((Pcal**2)/(Test**8.5))*(Lcal/(r**2))
        Deltaii1T = h*(fT-fT_values[81])
        Tcal = T_values[boundary]+h*fT-0.5*Deltaii1T
        
        # Executing Step 8
        ErelT = abs(Tcal-Test)/Tcal
        if ErelT < Erelmax:
            loop2 = False
        else:
            Test = Tcal
            Pest = Pcal
            
    # Executing Step 9
    nplus1_values[boundary+1] = (Tcal/Pcal)*((h*fP)/(h*fT))

    #####################################################################################################################################################

    # Boundary values (linear interpolation):
    r_boundary = []
    P_boundary = []
    T_boundary = []
    L_boundary = []
    M_boundary = []

    r_bound = ((r_values[boundary]-r_values[boundary+1])/(nplus1_values[boundary]-nplus1_values[boundary+1]))*(2.5-nplus1_values[boundary+1])+r_values[boundary+1]
    P_bound = ((P_values[boundary]-P_values[boundary+1])/(nplus1_values[boundary]-nplus1_values[boundary+1]))*(2.5-nplus1_values[boundary+1])+P_values[boundary+1]
    T_bound = ((T_values[boundary]-T_values[boundary+1])/(nplus1_values[boundary]-nplus1_values[boundary+1]))*(2.5-nplus1_values[boundary+1])+T_values[boundary+1]
    L_bound = ((L_values[boundary]-L_values[boundary+1])/(nplus1_values[boundary]-nplus1_values[boundary+1]))*(2.5-nplus1_values[boundary+1])+L_values[boundary+1]
    M_bound = ((M_values[boundary]-M_values[boundary+1])/(nplus1_values[boundary]-nplus1_values[boundary+1]))*(2.5-nplus1_values[boundary+1])+M_values[boundary+1]

    r_boundary.append(r_bound)
    P_boundary.append(P_bound)
    T_boundary.append(T_bound)
    L_boundary.append(L_bound)
    M_boundary.append(M_bound)
    
    #####################################################################################################################################################

    values_from_surface = pd.DataFrame ({
        'E': E_prod,
        'Phase': model_phase,
        'i': i_values,
        'r': r_values,
        'P': P_values,
        'T': T_values,
        'L': L_values,
        'M': M_values,
        'n+1': nplus1_values,
    })

    boundary_from_surface = pd.DataFrame ({
        'i boundary': boundary,
        'r boundary': r_boundary,
        'P boundary': P_boundary,
        'T boundary': T_boundary,
        'L boundary': L_boundary,
        'M boundary': M_boundary,   
    })

    return (K, values_from_surface, boundary_from_surface)
    