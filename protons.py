#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 10:09:17 2018

@author: mollykaplan

Calculates the proton spectrum, including ionization, pion and wind
advection losses. 

Outputs a plot of the spectrum versus energy. Also can output a tuple 
with the array of energy values in the first index and an array of 
the values of the proton spectrum at those energies in the second 
index. (This is commented out for now; using pSpec.py instead)
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as fun

#b for proton spectrum, including ionization, pion, and wind 
#advection losses. E is the kinetic energy.
def b(E):
    density=1e2
    beta=np.sqrt(1.-(E/.9383+1.)**(-2.))
    #return  (10+np.log(E)+E**(2.))
    if E>0.28172800:
        #includes pion decay since above threshhold
        return 5.775e-4*density*(10.9 + 2.*np.log(E/0.9383 + 1.)+ \
               2.*np.log(beta)-beta*beta)/beta+ 1.846e-2*density*E+E
    else:
        #pion losses not included
        return 5.775e-4*density*(10.9 + 2.*np.log(E/0.9383 + 1.)+ \
               2.*np.log(beta)-beta*beta)/beta+E


#inverse of the tau function. tau_0 and tau_c hard coded for now
def inv_tau(E):
    tau_0=1e-1
    return tau_0*np.sqrt(1.-(E/.9383+1.)**(-2.))*np.sqrt(E)+1.

#calculates integrand of integral in green's function
def func(E):
    return 1./b(E)*inv_tau(E)
 

#############################################################################
##Testing the number spectrum
Emin=1e-4
Emax=1e6
npts=100
x0=np.geomspace(Emin,Emax,npts) #equally spaced in log space
x=np.delete(np.delete(x0,0),-1) #delete first and last point

#array of number spectrum values
vec_Nspec=[]
for elt in x: vec_Nspec.append(fun.Nspec(func,b,elt,Emax,1.,-2.2))

#plotting the number spectrum
f1 = plt.figure(figsize=(11,6))
ax1 = f1.add_subplot(111)
plt.plot(x,vec_Nspec,Emin,Emax)
ax1.set_xscale('log')
ax1.set_yscale('log')
    
print(vec_Nspec)


def proton_spec():
    return (x,vec_Nspec)
