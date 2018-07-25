#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 10:47:02 2018

@author: mollykaplan

Calculates the injection spectrum of secondary electrons from an input
proton spectrum

Function secondary_e_spec outputs a tuple with the array of energy values 
in the first index and an array of the values of the secondary electron
spectrum at those energies in the second index. Takes energy values,
density and the timescale as inputs.

In testing section, outputs a plot of the spectrum and prints the values 
of the spectrum.
"""

import numpy as np 
import matplotlib.pyplot as plt
import functions as fun
from scipy.interpolate import interp1d

#important constants
mp=.9383   #proton mass [GeV/c^2]
me=5.1099891e-4   #electron mass [GeV/c^2]
s=me/mp   # m_e/(A_i*m_p)
Epmax=1e6

"""
Function that finds the knock-on electron cross section. Ee is
electron energy, and Ep is proton energy. Slightly modified version
of Torres 2004 Eq. 4.
"""
def ko_xsec(Ee,Ep):
    gamma_e=Ee/5.1099891e-4 + 1.
    gamma_p=Ep/0.9383 + 1.
    return 1./np.sqrt(1.-gamma_p**(-2.)) \
           *( (gamma_e - 1.)**(-2.)\
           -s*(gamma_p+.5*(s + 1./s))/((gamma_e - 1.)*gamma_p**(2.))\
           +.5*s*s/(gamma_p**(2.)))
           
"""
Function that returns the injection spectrum for secondary 
electrons at some Ee, as given by Eq. 5 from Torres 2004. 
interp_func is an interpolation over the proton spectrum. 
"""
def ko_source(Ee,density,interp_func):
    gamma_e=Ee/5.1099891e-4 + 1.
    E1p=mp*(.5*s*(gamma_e-1.)+\
        np.sqrt(1.+.5*(1.+s*s)*(gamma_e-1.)+\
        .25*s*s*(gamma_e-1.)**(2.))-1.)

    integrand= lambda Ep: ko_xsec(Ee,Ep)*10.**interp_func(np.log10(Ep))
    
    return 1616.5390*density*fun.fint(integrand,E1p,Epmax)

"""
Outputs a tuple with frequency values given by an array "x", and
secondary electron spectrum values at the values in x.
"""
def ko_sec(x,vec_p,NEp,density):
    #creating interpolated proton spectrum
    l_p=[]
    for elt in vec_p: l_p.append(np.log10(elt))
    l_NEp=[]
    for elt in NEp: l_NEp.append(np.log10(elt))
    interp_func=interp1d(l_p,l_NEp, 
                         fill_value=-np.Inf, 
                         bounds_error=False, 
                         kind='cubic')
    
    #creating array to hold injection spectrum data
    vec_KO=np.zeros(len(x))
    for i in range(len(x)): vec_KO[i]=ko_source(x[i],density,interp_func)
        
    return (x,vec_KO)


#Testing the secondary electron spectrum
##############################################################################
##creating energy array with points equally spaced in log space and
##deleting first and last point
#Emin=1e-4
#Emax=1e5
#npts=100
#x=np.delete(np.delete(np.geomspace(Emin,Emax,npts),0),-1) 
#
#density=1e3
#tau_0=1e-1
#
##creating array of injection spectrum values
#spec=secondary_e_spec(x,density,tau_0)
#vec_KO=spec[1]
#vec_x=spec[0]
#        
##Printing energy and knock-on spectrum values
#print("x vals: "+str(vec_x))
#print("K0 vals: "+str(vec_KO))
#    
#
##plotting the injection spectrum
#f1 = plt.figure(figsize=(11,6))
#ax1 = f1.add_subplot(111)
#plt.plot(vec_x,vec_KO)
#ax1.set_xscale('log')
#ax1.set_yscale('log')
