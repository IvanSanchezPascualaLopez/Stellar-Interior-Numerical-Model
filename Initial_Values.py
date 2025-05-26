# -*- coding: utf-8 -*-
"""
Created on Tue May  6 19:38:34 2025

@author: sanch
"""

# Necessary libraries:
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

#This functions asks for the initual values of the star:
def Initial_Values ():
    X = float(input("Insert a value for X (mass fraction of hydrogen) [adimensional]:"))
    Y = float(input("Insert a value for Y (mass fraction of helium) [adimensional]:"))
    M_tot = float(input("Insert a value for M_tot (total mass of the star) [10^33 g]:"))
    R_tot = float(input("Insert a value for R_tot (radius of the star) [10^10 cm]:"))
    L_tot = float(input("Insert a value for L_tot (luminosity of the star) [10^33 erg s^-1]:"))
    Tc = float(input("Insert a value for Tc (central temperature of the star) [10^7 K]:"))
    return X, Y, M_tot, R_tot, L_tot, Tc