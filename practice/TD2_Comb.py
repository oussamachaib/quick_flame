#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 22:23:29 2021

@author: oussamachaib
"""

from pylab import*
import matplotlib as plt
import pandas as pd

# molar fractions
MH=1
MC=12
MN=14
MO=16


# hydrocarbon structure (format: CnHm)
# usually, m = 2n+2
n=1
m=4
if(n==1):
    n2=""
else:
    n2=n


# air mixture (default: 3.76, pure oxy: 0)
beta=3.76
p=n+m/4


# reactants or products (1: reac, 3: prod)
rp=3

phi_range=arange(0,5.1,.1)
Y_fuel=zeros(len(phi_range))
Y_O2=zeros(len(phi_range))
Y_N2=zeros(len(phi_range))
Y_CO2=zeros(len(phi_range))
Y_H2O=zeros(len(phi_range))

i=0
for phi in phi_range:
    # mass fractions
    column_names=["nk_reac [ ]","Yk_reac [ ]","nk_prod [ ]","Yk_prod [ ]","M [g/mol]"]
    index_names=[f"C{n2}H{m}","O2","N2","CO2","H2O"]
    mf=pd.DataFrame(data=zeros((len(index_names),len(column_names))),index=index_names,columns=column_names)
    
    data=zeros((len(index_names),len(column_names)))
    mf["M [g/mol]"]=[(n*MC+m*MH),2*MO,2*MN,MC+2*MO,2*MH+MO]
    
    # reactants
    
    mf.iloc[:,0]=[phi,p,p*beta,0,0]
    if(phi<=1):
        mf.iloc[:,2]=[0,p*(1-phi),p*beta,phi*n,phi*m/2]
    else:
        mf.iloc[:,2]=[phi-1,0,p*beta,n,m/2]    
    
    # molar mass
    M_tot=sum(mf.iloc[:,-1]*mf.iloc[:,0])
    
    # calculating mass fractions
    mf.iloc[:,1]=mf.iloc[:,0]*mf.iloc[:,-1]/M_tot
    mf.iloc[:,3]=mf.iloc[:,2]*mf.iloc[:,-1]/M_tot

    Y_fuel[i]=mf.iloc[0,rp]
    Y_O2[i]=mf.iloc[1,rp]
    Y_N2[i]=mf.iloc[2,rp]
    Y_CO2[i]=mf.iloc[3,rp]
    Y_H2O[i]=mf.iloc[4,rp]

    i=i+1

figure(1)

plot(phi_range,Y_fuel)
plot(phi_range,Y_O2)
plot(phi_range,Y_N2)
plot(phi_range,Y_CO2)
plot(phi_range,Y_H2O)

legend(["Fuel",r"$O_2$",r"$N_2$",r"$CO_2$",r"$H_2O$"],bbox_to_anchor=(1,1), loc="upper left")
title(r'$Y_k = f(\phi)$')
show()
xlabel(r'$\phi$')
ylabel(r'$Y_k$')
tight_layout()





