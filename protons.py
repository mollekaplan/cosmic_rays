#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 10:09:17 2018

@author: mollykaplan

Calculates the proton spectrum, including ionization, pion and wind
advection losses.

Function proton_spec outputs a tuple with the array of energy values 
in the first index and an array of the values of the proton spectrum at 
those energies in the second index. Takes density and the timescale as inputs.
"""

import numpy as np
import functions as fun

"""
Energy loss rate for proton spectrum as given by Paglione et. al. 
1996 Eq. 7. Includes ionization, pion, and wind advection losses. 
E is the proton kinetic energy.
"""
def b(E,density):
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

"""
Inverse of the tau function as given by Paglione et. al. 2012, 
Eq. 3. The input parameter tau_0 is the inverse timescale.
"""
def inv_tau(E,tau_0):
    return tau_0*np.sqrt(1.-(E/.9383+1.)**(-2.))*np.sqrt(E)+1.

"""
Calculates integrand of integral in Green's function.
"""
def func(E,density,tau_0):
    return 1./b(E,density)*inv_tau(E,tau_0)



##############################################################################
"""
Outputs proton spectrum value array. Takes as input certain 
parameters of the ISM, including the ambient gas density and
the escape timescale.
"""
def proton_spec(density,tau_0):
    #proton energies for creating the proton spectrum array
    Emin=1e-4
    Emax=1e6
    npts=100
    #equally spaced in log space, deleting first and last point
    x=np.delete(np.delete(np.geomspace(Emin,Emax,npts),0),-1)
    
    #creating a proton spectrum array
    pspec=[]
    for elt in x: pspec.append(fun.Nspec_p(lambda E: func(E,density,tau_0),\
              lambda E: b(E,density),elt,Emax,1.,2.2))
    return (x,pspec)
