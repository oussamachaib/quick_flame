#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 20:25:05 2021

@author: oussamachaib
"""

from pylab import*
from matplotlib import*

# molar masses
WC=12
WH=1
WO=16
WN=14

# fuel mass or mass flow rate
mfuel=4
# air mass or mass flow rate 
mox=75

# hydrocarbon structure 
n=4
m=10

# mixture fraction at stoichiometry
alphas=(n*WC+m*WH)/((13/2)*(2*WO+3.76*2*WN))

# current mixture fraction
alpha=mfuel/mox

# equivalence ratio
phi=alpha/alphas

print(f'Equivalence ratio: {phi:.2}')