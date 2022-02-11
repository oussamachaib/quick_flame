#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 00:00:52 2021

@author: oussamachaib
"""

from pylab import*
import matplotlib.pyplot as plt

# hydrocarbon structure (format: CnHm)
# usually, m = 2n+2
n=8
m=18
if(n==1):
    n2=""
else:
    n2=n

# air mixture (default: 3.76, pure oxy: 0)
beta=3.76
p=n+m/4

# balancing equation
fig=plt.figure(figsize=(10,1))
plt.axis('off')
if(beta!=0):
    plt.text(0, 0.5,
         rf'Balanced equation:   $C_{{{n2}}}H_{{{m}}} + {p:.3}\,(O_2+{beta}\,N_2)$'
         +u' \u27F6 '
         +rf'${n}\,CO_2 + {m/2:.2}\,H_2O + {beta}\,({p})\,N_2$'
         ,fontsize=14)
else:
    plt.text(0, 0.5,
         rf'Balanced equation:   $C_{{{n}}}H_{{{m}}} + {p:.3}\,(O_2)$'
         +u' \u27F6 '
         +rf'${n}\,CO_2 + {m/2:.2}\,H_2O$'
         ,fontsize=14)

plt.show()
plt.tight_layout()
