# -*- coding: utf-8 -*-
"""
Created on Thu May 22 19:54:46 2025

@author: sanch
"""
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.font_manager as font_manager
import os

def Hipparcos_HR_Diagram (T_eff, L_tot):
    
    # The next lines are based on the code by Silva Lima F.J.: www.kaggle.com/code/fernandolima23/hrdiagram
    for dirname, _, filenames in os.walk('/kaggle/input'):
        for filename in filenames:
            print(os.path.join(dirname, filename))
    Data = pd.read_csv("hipparcos-voidmain.csv")
    pd.set_option("display.max_columns", 78)
    Data.head()
    
    New_dataframe = Data.loc[:,["Vmag", "Plx", "B-V"]][Data.loc[:,["Vmag", "Plx", "B-V"]]["Plx"]>0]
    New_dataframe.shape
    
    
    percentual_of_data_missing_data_missing = 100*(New_dataframe.isnull().sum()/len(New_dataframe["Plx"]))
    
    
    New_dataframe = New_dataframe.dropna(subset = ["B-V"])
    
    New_dataframe = New_dataframe.sample(frac=0.25) 
    
    New_dataframe.dtypes
    
    Magnitude_B_V = np.array(New_dataframe["B-V"], dtype=float)
    Magnitude_B_V2 = Magnitude_B_V.tolist()
    apparent_magnitude = np.array(New_dataframe["Vmag"], dtype=float)
    geometric_parallax = np.array(New_dataframe["Plx"], dtype=float)
    
    dispc = 1/geometric_parallax*1000
    type(dispc)
    
    absolute_magnitude = apparent_magnitude-5*np.log10(dispc/10)
    absolute_magnitude2 = absolute_magnitude.tolist()
    
    
    Font1 = {"family": "serif", "weight": "bold", "color": "black", "size": 13}
    Font2 = {"family": "serif", "weight": "bold", "color": "black", "size": 15}
    Font3 = font_manager.FontProperties(family="serif",
                                       weight='bold',
                                       style='normal', size=14) # Específico para a legenda
    
    
    Teff_values=4600*((1/(0.92*Magnitude_B_V+1.7))+(1/(0.92*Magnitude_B_V+0.62)))
    L_values=10**(-0.4*(absolute_magnitude-4.83))
    
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(Teff_values, np.log(L_values), c=Magnitude_B_V,  cmap = "autumn", marker = '.', s = 1, linewidth = 0.01)
    plt.plot(T_eff, np.log((L_tot*(10**33)*(10**-7))/3.8e26), marker='x', color='black', markersize=5, mew=1)
    plt.plot(5778, np.log(1), marker='x', color='black', markersize=5, mew=1)
    plt.plot(9662, np.log((75.904*(10**33)*(10**-7))/3.8e26), marker='x', color='black', markersize=5, mew=1)

    
    plt.text(T_eff-300, np.log((L_tot*(10**33)*(10**-7))/3.8e26) +0.1, 'Model Star 2', color='black', fontsize=11)
    plt.text(5778-100, np.log(1)+0.5, 'The Sun', color='black', fontsize=11)
    plt.text(9662-100, np.log((75.904*(10**33)*(10**-7))/3.8e26)+0.6, 'Model Star 1', color='black', fontsize=11)
    plt.xlabel("$T_{eff} [K]$", fontdict = Font1)
    plt.ylabel("$log(L / L_\odot)$", fontdict = Font1)
    plt.title("Hertzsprung-Russell Diagram", fontdict = Font2) 
    plt.xlim(0,15000)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis = "both", direction = "in", labelcolor='black', labelsize='xx-large', top = True, right = True)
    ax.tick_params(which='major', direction = "in", color='black', length=5, width = 1)
    ax.tick_params(which='minor', direction = "in", length=2, color='black', width = 0.7, top = True, right = True)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_color("black")
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1)
    Color_BC = plt.gca()
    Color_BC.set_facecolor("white")
    Color_BC.patch.set_alpha(1)
    plt.gca().invert_xaxis()
    plt.legend(frameon = False, prop = Font3, labelcolor = "white")
    plt.tight_layout
    
    
    plt.grid(True, color='gray', linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.show()
    
    plt.show()