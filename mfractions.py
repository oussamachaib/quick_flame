#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 13:06:15 2021

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
n=8
m=18
if(n==1):
    n2=""
else:
    n2=n


# air mixture (default: 3.76, pure oxy: 0)
beta=3.76
p=n+m/4

# equivalence ratio
phi=0.6

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

# molar fractions
molf=mf.copy()
molf.columns=["nk_reac [ ]","Xk_reac [ ]","nk_prod [ ]","Xk_prod [ ]","M [g/mol]"]
molf.index=[f"C{n2}H{m}","O2","N2","CO2","H2O"]

# calculating molar fractions
molf.iloc[:,1]=molf.iloc[:,0]/sum(molf.iloc[:,0])
molf.iloc[:,3]=molf.iloc[:,2]/sum(molf.iloc[:,2])

print('Mass fraction table: \n')
print(mf)
print('\n')
print('Molar fraction table: \n')
print(molf)
print('\n')



