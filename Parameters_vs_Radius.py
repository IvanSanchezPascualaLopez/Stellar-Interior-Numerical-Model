# -*- coding: utf-8 -*-
"""
Created on Thu May 22 19:26:17 2025

@author: sanch
"""
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from Energy_Production_Rate import Energy_Production_Rate

def Parameters_vs_Radius (X, Y, M_tot, R_tot, L_tot, Tc, star_values, boundary_values):

    # Star Parameters:
    ##########################################################################    
    plt.figure(figsize = (5, 5))
    
    # Relative Parameters:
    plt.plot(star_values['r']/R_tot, star_values['P']/max(star_values['P']), '-', color = 'blue', label = 'Pressure')
    plt.plot(star_values['r']/R_tot, star_values['T']/max(star_values['T']), '-', color = 'red', label = 'Temperature')
    plt.plot(star_values['r']/R_tot, star_values['L']/max(star_values['L']), '-', color = 'gold', label = 'Luminosity')
    plt.plot(star_values['r']/R_tot, star_values['M']/max(star_values['M']), '-', color = 'purple', label = 'Mass')
    
    # Title & labels
    plt.title(r'Star Parameters', fontsize = 15)
    plt.xlabel(r'$r/R_{tot}$', fontsize = 15)
    plt.ylabel(r'$Relative \ value$', fontsize = 15)
    
    # X Axix interval:
    plt.xlim(0, 1)
    
    # Tick axis
    plt.tick_params(axis = 'both', which = 'major', labelsize = 11)
    plt.tick_params(bottom = True, top = True, left = True, right = True)
    
    # Bordes de la figura
    plt.gca().spines[['bottom', 'left', 'top', 'right']].set_linewidth(1)
    
    # Regions:
    plt.axvspan(0, star_values.at[109, 'r']/R_tot, label = 'Center', color = 'indigo', alpha = 0.2, linestyle = '-')  # First three layers (convective)
    plt.axvspan(star_values.at[109, 'r']/R_tot, boundary_values.at[0, 'r boundary']/R_tot, label = 'Convective Core', color = 'royalblue', alpha = 0.2, linestyle = '-')  # Convective Core
    plt.axvspan(boundary_values.at[0, 'r boundary']/R_tot, star_values.at[51, 'r']/R_tot, label = 'Radiative Envelope', color = 'aquamarine', alpha=0.1, linestyle = '-')  # Radiative Envelope
    plt.axvspan(star_values.at[51, 'r']/R_tot, star_values.at[14, 'r']/R_tot, label = 'pp cycle starts', color = 'palegreen', alpha = 0.2, linestyle = '-')  # pp cycle starts
    plt.axvspan(star_values.at[14, 'r']/R_tot, star_values.at[11, 'r']/R_tot, label = 'Beginning', color = 'yellow', alpha = 0.2, linestyle = '-')  # First three layers (radiative)
    plt.axvspan(star_values.at[11, 'r']/R_tot, 1, label = 'Outer Layers', color = 'coral', alpha = 0.2, linestyle = '-')  # Outer layers
    
    # Legend:
    plt.legend(fontsize = 11, markerscale = 1)
    
    # Grid
    plt.grid(True, alpha = 0.1)
    
    plt.show()
    ###########################################################################
    #Density:
    #########################################################################
    Z = 1-X-Y                    # Mass fraction of heavy elements                           [adimensional]
    mu = (1)/(2*X+3*Y/4+Z/2)     # Average molecular weight (constant throughout the star)   [adimensional]
    Na = 6.022e23                # Avogadro number                                           [adimensional]
    H = 1/Na                     # Avogadro numer inverse                                    [adimensional]
    k = 1.380649e-23             # Boltzmann constant                                        [din*cm/10^7K]=[J/K]
    
    #Density values:
    rho = ((mu*H*star_values['P'])/(k*star_values['T']))       # Density 
    print('The central density of the star, rho_c, is', round(max(rho)*100, 3), 'g/cm^3')
    
    
    # Plot:
    plt.figure(figsize = (5, 5))
    
    
    plt.plot(star_values['r']/R_tot, rho/max(rho), '-', color = 'green', label = 'Density')
    
    # Title & labels:
    plt.title(r'Density', fontsize = 15)
    plt.xlabel(r'$r/R_{tot}$', fontsize = 15)
    plt.ylabel(r'$\rho/\rho_c$', fontsize = 15)
    
    # X axis interval:
    plt.xlim(0, 1)
    
    # Tick axis:
    plt.tick_params(axis = 'both', which = 'major', labelsize = 11)
    plt.tick_params(bottom = True, top = True, left = True, right = True)
    
    # Regions:
    plt.axvspan(0, star_values.at[109, 'r']/R_tot, label = 'Center', color = 'indigo', alpha = 0.2, linestyle = '-')  # First three layers (convective)
    plt.axvspan(star_values.at[109, 'r']/R_tot, boundary_values.at[0, 'r boundary']/R_tot, label = 'Convective Core', color = 'royalblue', alpha = 0.2, linestyle = '-')  # Convective Core
    plt.axvspan(boundary_values.at[0, 'r boundary']/R_tot, star_values.at[51, 'r']/R_tot, label = 'Radiative Envelope', color = 'aquamarine', alpha=0.1, linestyle = '-')  # Radiative Envelope
    plt.axvspan(star_values.at[51, 'r']/R_tot, star_values.at[14, 'r']/R_tot, label = 'pp cycle starts', color = 'palegreen', alpha = 0.2, linestyle = '-')  # pp cycle starts
    plt.axvspan(star_values.at[14, 'r']/R_tot, star_values.at[11, 'r']/R_tot, label = 'Beginning', color = 'yellow', alpha = 0.2, linestyle = '-')  # First three layers (radiative)
    plt.axvspan(star_values.at[11, 'r']/R_tot, 1, label = 'Outer Layers', color = 'coral', alpha = 0.2, linestyle = '-')  # Outer layers
    
    # Legend:
    plt.legend(fontsize = 11, markerscale = 1)
    
    # Grid:
    plt.grid(True, alpha = 0.1)
    
    plt.show()
    ###########################################################################
    #Energy Production rate:
    #######################################################################
    r = star_values['r']
    T = star_values['T']
    P = star_values['P']
    eps1 = np.zeros(len(r))
    nu = np.zeros(len(r))
    X1 = np.zeros(len(r))
    X2 = np.zeros(len(r))
    eps_values = []
    for i in range (0, len(r)):
        cycle1, eps1[i], nu[i], X1[i], X2[i] = Energy_Production_Rate (T[i], P[i], X, Y)
        eps_values.append(eps1[i]*X1[i]*X2[i]*rho[i]*100*((T[i]*10)**nu[i]))
    
    # Plot:
    plt.figure(figsize = (5, 5))
    
    plt.plot(r/R_tot, eps_values/max(eps_values), '-', color='orange', label = 'Energy Production Rate')
    
    # Title & axis:
    plt.title(r'Energy production rate', fontsize = 15)
    plt.xlabel(r'$r/R_{tot}$', fontsize = 15)
    plt.ylabel(r'$\epsilon_{rel}$', fontsize = 15)
    
    # X interval:
    plt.xlim(0, 1)
    
    # Tick parameters:
    plt.tick_params(axis = 'both', which = 'major', labelsize = 11)
    plt.tick_params(bottom = True, top = True, left = True, right = True)
    
    # Regions:
    plt.axvspan(0, star_values.at[109, 'r']/R_tot, label = 'Center', color = 'indigo', alpha = 0.2, linestyle = '-')  # First three layers (convective)
    plt.axvspan(star_values.at[109, 'r']/R_tot, boundary_values.at[0, 'r boundary']/R_tot, label = 'Convective Core', color = 'royalblue', alpha = 0.2, linestyle = '-')  # Convective Core
    plt.axvspan(boundary_values.at[0, 'r boundary']/R_tot, star_values.at[51, 'r']/R_tot, label = 'Radiative Envelope', color = 'aquamarine', alpha=0.1, linestyle = '-')  # Radiative Envelope
    plt.axvspan(star_values.at[51, 'r']/R_tot, star_values.at[14, 'r']/R_tot, label = 'pp cycle starts', color = 'palegreen', alpha = 0.2, linestyle = '-')  # pp cycle starts
    plt.axvspan(star_values.at[14, 'r']/R_tot, star_values.at[11, 'r']/R_tot, label = 'Beginning', color = 'yellow', alpha = 0.2, linestyle = '-')  # First three layers (radiative)
    plt.axvspan(star_values.at[11, 'r']/R_tot, 1, label = 'Outer Layers', color = 'coral', alpha = 0.2, linestyle = '-')  # Outer layers
    
    # Legend:
    plt.legend(fontsize = 11, markerscale = 1)
    
    # Grid:
    plt.grid(True, alpha = 0.1)
    
    plt.show()
    ###########################################################################
    #Opacity:
    #########################################################################
    kappa = np.zeros(len(r))
    for i in range(1, len(r)):
        kappa[i] = 4.34*1e25*Z*(X+1)*rho[i]/(T[i]**3.5)
    
    # Plot:
    plt.figure(figsize=(5, 5))
    
    plt.plot(r/R_tot, kappa/max(kappa), '-', color = 'brown', label = 'Opacity')
    
    # Title & labels:
    plt.title(r'Opacity', fontsize = 15)
    plt.xlabel(r'$r/R_{tot}$', fontsize = 15)
    plt.ylabel(r'$\kappa_{rel}$', fontsize = 15)
    
    # X, Y intervals:
    plt.xlim(0, 1)
    plt.ylim(0, 1.05)
    
    # Tick parameters
    plt.tick_params(axis = 'both', which = 'major', labelsize = 11)
    plt.tick_params(bottom = True, top = True, left = True, right = True)
    
    # Regions:
    plt.axvspan(0, star_values.at[109, 'r']/R_tot, label = 'Center', color = 'indigo', alpha = 0.2, linestyle = '-')  # First three layers (convective)
    plt.axvspan(star_values.at[109, 'r']/R_tot, boundary_values.at[0, 'r boundary']/R_tot, label = 'Convective Core', color = 'royalblue', alpha = 0.2, linestyle = '-')  # Convective Core
    plt.axvspan(boundary_values.at[0, 'r boundary']/R_tot, star_values.at[51, 'r']/R_tot, label = 'Radiative Envelope', color = 'aquamarine', alpha=0.1, linestyle = '-')  # Radiative Envelope
    plt.axvspan(star_values.at[51, 'r']/R_tot, star_values.at[14, 'r']/R_tot, label = 'pp cycle starts', color = 'palegreen', alpha = 0.2, linestyle = '-')  # pp cycle starts
    plt.axvspan(star_values.at[14, 'r']/R_tot, star_values.at[11, 'r']/R_tot, label = 'Beginning', color = 'yellow', alpha = 0.2, linestyle = '-')  # First three layers (radiative)
    plt.axvspan(star_values.at[11, 'r']/R_tot, 1, label = 'Outer Layers', color = 'coral', alpha = 0.2, linestyle = '-')  # Outer layers
    
    # Legend:
    plt.legend(fontsize = 11, markerscale = 1)
    
    # Grid
    plt.grid(True, alpha=0.1)
    
    plt.show()
    #############################################################################
    
    
    return rho, eps_values, kappa