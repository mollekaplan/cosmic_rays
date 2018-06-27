#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 13:26:44 2018

@author: mollykaplan

This code calculates the synchrotron emission from an input electron
spectrum. Outputs plots of the emissivity versus frequency (nu) 
for 20 (nfreq) different values of nu.
"""

import numpy as np
from e_spec import electron_spec
import matplotlib.pyplot as plt
import functions as fun
from bessel import Fx
#from bessel import mod_bessel 

vec_e=electron_spec()[0]  #array of energy values used inelectron spectrum
NEe=electron_spec()[1]  #array of electron spectrum values

x=Fx()[0]  #array of x values used for bessel integral
F=Fx()[1]  #array of bessel integral values

mag=1e3 #B-field in microGauss


#function which calulates the integrand of the energy integration, which 
#includes the electron spectrum, the bessel integral, and the integral 
#over pitch angle (theta)
def integrand(E,nu):
    Nspec_e=fun.log_interp(vec_e,NEe,E) #electron spectrum at energy E
    
    x0=0.23815917*nu/mag  #nu/nu_c * gamma^2 * sin(theta)
    g2=(E/5.1099891e-4 + 1.)**2.  #gamma^2
    min_th=np.arcsin(x0*.01/g2)  #ensures nu/nu_c <= 100
    
    dth=(np.pi/2. - min_th)*.1  #division width 
    th=np.linspace(min_th,np.pi/2.,11)  #array of 11 theta values
  
    #calculating sin^2(theta) * F(xj) at each theta value
    sinF=[]
    for theta in th:
        xj=x0/(g2*np.sin(theta))  #nu/nu_c
        sinF.append(np.sin(theta)**2.*fun.log_interp(x,F,xj,True))
    
    #trapezoidal integration over theta
    thint=dth*(.5*(sinF[0]+sinF[10]))
    for i in range(1,10): thint+=dth*sinF[i]
    
    return 1.86558e-6*mag*Nspec_e*thint

    
def emissivity(nu):    
    x0=0.23815917*nu/mag
    minE=5.1099891e-4*(np.sqrt(.01*x0)-1.) #ensures sin(theta)<=1 for nu/nu_c<=100

    #values where the integral is split up
    endpts=[1e-2,3e-2,4e-2,5e-2,6e-2,6.5e-2,7.5e-2,8.5e-2,1e-1,\
            1.2e-1,1.3e-1,1.5e-1,1.7e-1,1.9e-1,2.1e-1,2.5e-1,3e-1,\
            3.2e-1,3.5e-1,4e-1,4.5e-1,5e-1,5.5e-1,6.5e-1,7.5e-1,\
            8.5e-1,1e0,1.3e0,1.7e0,2.2e0,3e0,5e0,8e0,1e1,1e2] 
     
    #integration over E
    return fun.segmented_int(lambda E: integrand(E,nu),\
            endpts,minE,Emax)

##############################################################################
#creating array of frequency values
Emin=1e-4
Emax=1e6
numin=1e8
numax=10**(11.18)
nfreq=20
nu=np.geomspace(numin,numax,nfreq)

#creating vector of emissivity values and printing out nu and ev(nu)
ev=[]
for elt in nu: 
    print(str(elt)+"--"+str(emissivity(elt)))
    ev.append(emissivity(elt))

##plotting the synchrotron spectrum
f1 = plt.figure(figsize=(11,6))
ax1 = f1.add_subplot(111)
plt.plot(nu,ev(nu),numin,numax)
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.ylim(2e-7,6.3e-5)

#printing out values of emissivity and minimum energy at frequency nu_val
nu_val=1e8
print(emissivity(nu_val))
x0=0.23815917*nu_val/mag
minE=5.1099891e-4*(np.sqrt(.01*x0)-1.)
print(minE)



