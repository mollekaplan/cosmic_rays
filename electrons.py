#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 09:56:08 2018

@author: mollykaplan

Calculates the electron spectrum, including ionization, bremsstrahlung, 
wind advection, inverse Compton and synchrotron losses. Pulls the energy loss 
rate in the Klein-Nishina limit from the function b_ic in ic.py.

Function electron_spec outputs a tuple with the array of energy values 
in the first index and an array of the values of the proton spectrum at 
those energies in the second index. Takes density, magnetic field magnitude, 
energy density and the timescale as inputs.
"""

import numpy as np
import functions as fun
from scipy.interpolate import interp1d
from ic import b_ic
from secondary_e import secondary_e_spec

"""
The energy loss rate, b, for the electron spectrum; a slightly 
modified version of Paglione et. al. 1996 Eq. 2 to include 
inverse Compton losses in the Klein-Nishina limit. Includes 
ionization, bremsstrahlung, wind advection, inverse Compton and 
synchrotron losses. E is the kinetic energy. ic_interp is the 
function of IC losses; not included in b to reduce run-times.
"""     
def b(E,density,mag,Uph,ic_interp):
    mag2=mag*mag
    
    #find IC contribution to the energy loss rate
    ic_contr= 10.**ic_interp(np.log10(E))
    
    return 8.678e-4*density*(6.85 + np.log(E/5.1099891e-4 + 1.)+\
           38.305268*E) + E + 7.89e-5*mag2*(E + 5.1099891e-4)**2. + \
           ic_contr

"""
Inverse of the tau function as given by Paglione et. al. 2012, 
Eq. 3. The input parameter tau_0 is the inverse timescale.
"""
def inv_tau(E,tau_0):
    g2=(E/5.1099891e-4 + 1.)**2.
    return tau_0*np.sqrt(E - E/g2)+1.

"""
Calculates integrand of integral in Green's function.
"""
def func(E,density,mag,Uph,tau_0,ic_interp):
    return 1./b(E,density,mag,Uph,ic_interp)*inv_tau(E,tau_0)



##############################################################################
"""
Outputs electron spectrum value array. Takes as input certain 
parameters of the ISM, including the ambient gas density, 
the magnetic field magnitude, the photon energy density, and
the escape timescale.
"""
def electron_spec(density,mag,Uph,tau_0):
    #electron enrgies for creating the electron spectrum array
    Emin=1e-4
    Emax=1e6
    npts=100
    #equally spaced in log space, deleting first and last point
    x=np.delete(np.delete(np.geomspace(Emin,Emax,npts),0),-1)
    
    #creating electron energy spectrum for IC contributions
    Emin_ic=1e-4
    Emax_ic=1e6
    npts_ic=100
    Ee=np.geomspace(Emin_ic,Emax_ic,npts_ic)

    #finding IC contribution to energy loss rate
    ic=b_ic(Ee,Uph)
    l_E=[]
    for elt in Ee: l_E.append(np.log10(elt))
    l_ic=[]
    for elt in ic: l_ic.append(np.log10(elt))
    ic_interp=interp1d(l_E,l_ic, 
                  fill_value=-np.Inf, 
                  bounds_error=False, 
                  kind='cubic')
    
    #creating electron energy spectrum for secondary contributions
    Emin_s=1e-4
    Emax_s=1e5
    npts_s=100
    x_s=np.delete(np.delete(np.geomspace(Emin_s,Emax_s,npts_s),0),-1) 

    #finding the secondary lepton spectrum
    s_E=secondary_e_spec(x_s,density,tau_0)[0]
    s_spec=secondary_e_spec(x_s,density,tau_0)[1]
    l_s_E=[]
    for elt in s_E: l_s_E.append(np.log10(elt))
    l_s_spec=[]
    for elt in s_spec: l_s_spec.append(np.log10(elt))
    scndy_spec=interp1d(l_s_E,l_s_spec, 
                  fill_value=-np.Inf, 
                  bounds_error=False, 
                  kind='cubic')
    
    #creating an electron spectrum array
    espec=[]
    for elt in x: espec.append(fun.Nspec_e(lambda E: \
              func(E,density,mag,Uph,tau_0,ic_interp),
              lambda E: b(E,density,mag,Uph,ic_interp),
              elt,Emax,1.,2.2,scndy_spec))
    
    return (x,espec)



