#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 09:56:08 2018

@author: mollykaplan

Calculates the electron spectrum, including ionization, bremsstrahlung, 
wind advection, and synchrotron losses. Non-relativistic IC losses 
included, but Klein-Nishina limit is not yet considered.

Outputs a plot of the spectrum versus energy
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as fun

#b for electron spectrum, including ionization, bremsstrahlung, wind 
#advection, and synchrotron losses. Non-relativistic IC losses 
#are also included. E is the kinetic energy. 
def b(E):
    density=1e2
    mag2=1e4
    Uph=2e2
    return 8.678e-4*density*(6.85 + np.log(E/5.1099891e-4 + 1.)\
           + 38.305268*E) + E + (3.22e-3*Uph + 7.89e-5*mag2)*\
           (E + 5.1099891e-4)**2.

#inverse of the tau function. tau_0 and tau_c hard coded for now
def inv_tau(E):
    tau_0=1e-1
    g2=(E/5.1099891e-4 + 1.)**2.
    return tau_0*np.sqrt(E - E/g2)+1.

#calculates integrand of integral in green's function
def func(E):
    return 1./b(E)*inv_tau(E)

##############################################################################
##Testing the number spectrum
Emin=1e-4
Emax=1e6
npts=100
x0=np.geomspace(Emin,Emax,npts) #equally spaced in log space
x=np.delete(np.delete(x0,0),-1) #delete first and last point


#array of number spectrum values
vec_Nspec=[]
for elt in x: vec_Nspec.append(fun.Nspec(func,b,elt,Emax,1.,-2.2))

print(vec_Nspec) #print out spectrum

#plotting the number spectrum
f1 = plt.figure(figsize=(11,6))
ax1 = f1.add_subplot(111)
plt.plot(x,vec_Nspec,Emin,Emax)
ax1.set_xscale('log')
ax1.set_yscale('log')