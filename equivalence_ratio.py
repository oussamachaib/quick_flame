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
mfuel=.44
# oxidizer mass or mass flow rate (air or oxygen-enriched air)
mox=10.6

# hydrocarbon structure (format: CnHm)
n=1
m=4

# air mixture (default: 3.76, pure oxy: 0)
beta=0

# mixture fraction at stoichiometry
alphas=(n*WC+m*WH)/((n+m/4)*(2*WO+beta*2*WN))

# current mixture fraction
alpha=mfuel/mox

# equivalence ratio
phi=alpha/alphas

print(f'Equivalence ratio: {phi:.2}')