#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 22:23:29 2021

@author: oussamachaib
"""

from pylab import*
import matplotlib as plt
import pandas as pd
from scipy import interpolate
from sklearn.linear_model import LinearRegression

close('all')

#%% Constants

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


# air mixture
beta=3.76
p=n+m/4

# equivalence ratio
phi_range=arange(0,5.1,.1)

#%% Reactants mass fraction = f(phi)

rp=1

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

legend([fr"$C_{{{n2}}} H_{{{m}}}$",r"$O_2$",r"$N_2$"],bbox_to_anchor=(1,1), loc="upper left")
title(r'$Y_k = f(\phi)$')
show()
xlabel(r'$\phi$')
ylabel(r'$Y_k$')
suptitle('Reactants')
tight_layout()

#%% Products mass fraction = f(phi)
rp=3

Y_fuel=zeros(len(phi_range))
Y_O2=zeros(len(phi_range))
Y_N2=zeros(len(phi_range))
Y_CO2=zeros(len(phi_range))
Y_H2O=zeros(len(phi_range))

i=0
for phi in phi_range:
    # mass fractions
    column_names=["nk_reac [ ]","Yk_reac [ ]","nk_prod [ ]","Yk_prod [ ]","M [g/mol]"]
    index_names=[f"C_{n2}H_{m}","O2","N2","CO2","H2O"]
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

figure(2)

plot(phi_range,Y_fuel)
plot(phi_range,Y_O2)
plot(phi_range,Y_N2)
plot(phi_range,Y_CO2)
plot(phi_range,Y_H2O)

legend([fr"$C_{{{n2}}} H_{{{m}}}$",r"$O_2$",r"$N_2$",r"$CO_2$",r"$H_2O$"],bbox_to_anchor=(1,1), loc="upper left")
title(r'$Y_k = f(\phi)$')
show()
xlabel(r'$\phi$')
ylabel(r'$Y_k$')
suptitle('Products')
tight_layout()

#%% Adiabatic flame temperature = f(phi)

# Hyp 0: Cp,k = Cp,N2(T0)
Cp_N2=1039.285714 #J/kg
T0=298 
T2=[[] for i in range(4)]
i=0
for phi in phi_range:
    # mass fractions
    column_names=["nk_reac [ ]","Yk_reac [ ]","nk_prod [ ]","Yk_prod [ ]","M [g/mol]"]
    index_names=[f"C_{n2}H_{m}","O2","N2","CO2","H2O"]
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
    
    mf["dh0_f [J/kg]"]=[-4.7e06,0,0,-8.94e06,-1.35e07]
    mf["Cp(T0)"]=ones(5)*Cp_N2
    
    T2[0].append((sum(mf.iloc[:,1]*mf["dh0_f [J/kg]"])-sum(mf.iloc[:,3]*mf["dh0_f [J/kg]"]))/Cp_N2+T0)
    i=i+1

# Hyp 1: different species Cp but all at T0

T0=298
i=0
for phi in phi_range:
    # mass fractions
    column_names=["nk_reac [ ]","Yk_reac [ ]","nk_prod [ ]","Yk_prod [ ]","M [g/mol]"]
    index_names=[f"C_{n2}H_{m}","O2","N2","CO2","H2O"]
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
    
    mf["dh0_f [J/kg]"]=[-4.7e06,0,0,-8.94e06,-1.35e07]
    mf["Cp(T0) [J/kg.K]"]=[2237.5,915.625,1039.285714,843.1818182,1866.666667]
    
    T2[1].append((sum(mf.iloc[:,1]*mf["dh0_f [J/kg]"])-sum(mf.iloc[:,3]*mf["dh0_f [J/kg]"]))/(sum(mf.iloc[:,3]*mf["Cp(T0) [J/kg.K]"]))+T0)
    i=i+1

# Hyp 2: different species Cp but all at T=1000K

T0=298
i=0
for phi in phi_range:
    # mass fractions
    column_names=["nk_reac [ ]","Yk_reac [ ]","nk_prod [ ]","Yk_prod [ ]","M [g/mol]"]
    index_names=[f"C_{n2}H_{m}","O2","N2","CO2","H2O"]
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
    
    mf["dh0_f [J/kg]"]=[-4.7e06,0,0,-8.94e06,-1.35e07]
    mf["Cp(T0) [J/kg.K]"]=[4500,1087.5,1164.285714,1236.363636,2288.888889]
    
    T2[2].append((sum(mf.iloc[:,1]*mf["dh0_f [J/kg]"])-sum(mf.iloc[:,3]*mf["dh0_f [J/kg]"]))/(sum(mf.iloc[:,3]*mf["Cp(T0) [J/kg.K]"]))+T0)
    i=i+1

# Hyp 3: using tabulated sensible enthalpy

T0=298.15

temps=[298.15,600,2200,2300,2400,2500,2600,2700]
sh=pd.DataFrame(index=temps,columns=[f"C_{n2}H_{m}","O2","CO2","H2O","N2"])
sh.iloc[0,:]=zeros(5)
sh.iloc[1,:]=[13.1,9.2,12.9,10.5,9.9]
sh.iloc[2,:]=[142.7,66.8,103.5,83.1,63.4]
sh.iloc[3,:]=[152.4,70.6,109.6,88.4,67]
sh.iloc[4,:]=[162.1,74.4,115.8,93.7,70.6]
sh.iloc[5,:]=[172,78.3,122,99.1,74.3]
sh.iloc[6,:]=[181.9,82.3,128.1,104.5,78]
sh.iloc[7,:]=[191.9,86.1,134.2,110,81.6]

ffuel=interpolate.interp1d(temps, sh[f"C_{n2}H_{m}"])
fO2=interpolate.interp1d(temps, sh["O2"])
fCO2=interpolate.interp1d(temps, sh["CO2"])
fH2O=interpolate.interp1d(temps, sh["H2O"])
fN2=interpolate.interp1d(temps, sh["N2"])

i=0

for phi in phi_range:
    # mass fractions
    column_names=["nk_reac [ ]","Yk_reac [ ]","nk_prod [ ]","Yk_prod [ ]","M [g/mol]"]
    index_names=[f"C_{n2}H_{m}","O2","N2","CO2","H2O"]
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
    
    mf["dh0_f [J/kg]"]=[-4.7e06,0,0,-8.94e06,-1.35e07]
    mf["Cp(T0) [J/kg.K]"]=[2237.5,915.625,1039.285714,843.1818182,1866.666667]

    hsfuel=ones(len(temps))
    hsO2=ones(len(temps))
    hsCO2=ones(len(temps))
    hsH2O=ones(len(temps))
    hsN2=ones(len(temps))    
    for i in range(1,len(temps)+1):
        hsfuel[i-1]=trapz(ffuel(temps[0:i]),temps[0:i])
        hsO2[i-1]=trapz(fO2(temps[0:i]),temps[0:i])
        hsCO2[i-1]=trapz(fCO2(temps[0:i]),temps[0:i])
        hsH2O[i-1]=trapz(fH2O(temps[0:i]),temps[0:i])
        hsN2[i-1]=trapz(fN2(temps[0:i]),temps[0:i])
        
    rhs=pd.DataFrame(index=temps,columns=[f"C_{n2}H_{m}","O2","CO2","H2O","N2"])
    rhs[f"C_{n2}H_{m}"]=hsfuel*mf.iloc[0,4]*mf.iloc[0,3]
    rhs["O2"]=hsO2*mf.iloc[1,4]*mf.iloc[1,3]
    rhs["N2"]=hsN2*mf.iloc[2,4]*mf.iloc[2,3]
    rhs["CO2"]=hsCO2*mf.iloc[3,4]*mf.iloc[3,3]
    rhs["H2O"]=hsN2*mf.iloc[4,4]*mf.iloc[4,3]
    
    rhs["sum"]=rhs.sum(axis=1)
    
    frhs=interpolate.interp1d(rhs["sum"],temps)
    T2[3].append(frhs((sum(mf.iloc[:,1]*mf["dh0_f [J/kg]"])-sum(mf.iloc[:,3]*mf["dh0_f [J/kg]"]))))

figure(3)
for i in range(4):
    plot(phi_range,T2[i])
    
legend(["Hyp 0","Hyp 1","Hyp 2","Hyp 3"],bbox_to_anchor=(1,1), loc="upper left")
title(r'$T_2 = f(\phi)$')
show()
xlabel(r'$\phi$')
ylabel(r'$T_2$')
suptitle('Adiabatic flame temperature')
tight_layout()  