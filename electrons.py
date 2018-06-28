#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 09:56:08 2018

@author: mollykaplan

Calculates the electron spectrum, including ionization, bremsstrahlung, 
wind advection, and synchrotron losses. Non-relativistic IC losses 
included, but Klein-Nishina limit is not yet considered.

Function electron_spec outputs a tuple with the array of energy values 
in the first index and an array of the values of the proton spectrum at 
those energies in the second index. Takes density, magnetic field magnitude, 
energy density and the timescale as inputs.

Under testing section, can output a plot of the spectrum versus energy.
"""

import numpy as np
import matplotlib.pyplot as plt
import functions as fun

#b for electron spectrum, including ionization, bremsstrahlung, wind 
#advection, and synchrotron losses. Non-relativistic IC losses 
#are also included. E is the kinetic energy. 
def b(E,density,mag,Uph):
    mag2=mag*mag
    return 8.678e-4*density*(6.85 + np.log(E/5.1099891e-4 + 1.)\
           + 38.305268*E) + E + (3.22e-3*Uph + 7.89e-5*mag2)*\
           (E + 5.1099891e-4)**2.


#inverse of the tau function. tau_0 and tau_c hard coded for now
def inv_tau(E,tau_0):
    g2=(E/5.1099891e-4 + 1.)**2.
    return tau_0*np.sqrt(E - E/g2)+1.


#calculates integrand of integral in green's function
def func(E,density,mag,Uph,tau_0):
    return 1./b(E,density,mag,Uph)*inv_tau(E,tau_0)


##Creating energy values array and function which outputs spectrum
##############################################################################
Emin=1e-4
Emax=1e6
npts=100
#equally spaced in log space, deleting first and last point
x=np.delete(np.delete(np.geomspace(Emin,Emax,npts),0),-1)


#Outputs electron spectrum values using the given inputs
def electron_spec(density,mag,Uph,tau_0):
    espec=[]
    for elt in x: espec.append(fun.Nspec(lambda E: func(E,density,mag,Uph,tau_0),\
              lambda E: b(E,density,mag,Uph),elt,Emax,1.,-2.2))
    return (x,espec)


###Testing the number spectrum
###############################################################################
#density=1e3
#mag=1.5e2
#Uph=2.7e2
#tau_0=1e-1
#
#vec_Nspec=electron_spec(density,mag,Uph,tau_0)[1]
#print(vec_Nspec) #print out spectrum
#
##plotting the number spectrum
#f1 = plt.figure(figsize=(11,6))
#ax1 = f1.add_subplot(111)
#plt.plot(x,vec_Nspec,Emin,Emax)
#ax1.set_xscale('log')
#ax1.set_yscale('log')



