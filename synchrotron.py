#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 13:26:44 2018

@author: mollykaplan

Calculates the synchrotron emission from an input electron
spectrum. 

Function synch outputs an array of synchrotron emissivity values. 
Takes frequencies, density, magnetic field magnitude, energy density 
and the timescale as inputs.

NB: May not run since electrons.py has been updated since using
this code.
"""

import numpy as np
from electrons import electron_spec
import matplotlib.pyplot as plt
import functions as fun
from bessel import Fx 

x=Fx()[0]  #array of x values used for bessel integral
F=Fx()[1]  #array of bessel integral values

Emax=1e6

##############################################################################
"""
Function which calulates the integrand of the energy integration, 
as given by Marscher et. al. 1978 Eq. 15, which includes the 
electron spectrum, the bessel integral, and the integral over 
pitch angle (theta).
"""
def integrand(E,nu,mag,vec_e,NEe):
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

"""
Calculates the emissivity by integrating over energy values.
"""
def emissivity(nu,mag,vec_e,NEe):    
    x0=0.23815917*nu/mag
    minE=5.1099891e-4*(np.sqrt(.01*x0)-1.) #ensures sin(theta)<=1 for nu/nu_c<=100
    
    #values where the integral is split up
    endpts=[1e-2,1.5e-2,2e-2,3e-2,4e-2,5e-2,5.5e-2,6e-2,6.5e-2,7e-2,
            7.5e-2,8.5e-2,1e-1,1.2e-1,1.3e-1,1.5e-1,1.7e-1,1.9e-1,
            2.1e-1,2.5e-1,3e-1,3.2e-1,3.5e-1,4e-1,4.5e-1,5e-1,5.5e-1,
            6.5e-1,7.5e-1,8.5e-1,1e0,1.3e0,1.7e0,2.2e0,3e0,5e0,8e0,
            1e1,1e2] 
     
    #integration over E
    return fun.segmented_int(lambda E: \
            integrand(E,nu,mag,vec_e,NEe),\
            endpts,minE,Emax)

"""
Outputs synchrotron emissivity values at the frequency values 
given by the array "nu"
"""
def synch(nu,density,mag,Uph,tau_0,s):
    pts=electron_spec(density,mag,Uph,tau_0,s)
    vec_e=pts[0]  #array of energy values used in electron spectrum
    NEe=pts[1]  #array of electron spectrum values
    
    synch_ev=[]
    for elt in nu: synch_ev.append(emissivity(elt,mag,vec_e,NEe))
    return synch_ev

"""
Function for testing a power law electron spectrum.
"""
def synch_pl(x,nu,mag,pl):
    arr=[]
    for elt in x: arr.append(elt**(pl))
    
    synch_ev=[]
    for elt in nu: synch_ev.append(emissivity(elt,mag,x,arr))
    return synch_ev
