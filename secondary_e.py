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
from protons import proton_spec #change to proton when done testing
import matplotlib.pyplot as plt
import functions as fun

#important constants
re=2.82e-13  #electron radius [cm]
mp=.9383   #proton mass [GeV/c^2]
me=5.1099891e-4   #electron mass [GeV/c^2]
s=me/mp   # m_e/(A_i*m_p)
c=9.4542e23   #cm Myr^(-1)

Epmax=1e6

#function that finds the knock-on electron cross section. Ee is
#electron energy, and Ep s proton eneergy
def KO_xsec(Ee,Ep):
    gamma_e=Ee/5.1099891e-4 + 1.
    gamma_p=Ep/0.9383 + 1.
    return 1./np.sqrt(1.-gamma_p**(-2.)) \
           *( (gamma_e - 1.)**(-2.)\
           -s*(gamma_p+.5*(s + 1./s))/((gamma_e - 1.)*gamma_p**(2.))\
           +.5*s*s/(gamma_p**(2.)))
           

#function that returns the injection spectrum for secondary 
#electrons at some Ee.
def KO_source(Ee,density,vec_p,NEp):
    gamma_e=Ee/5.1099891e-4 + 1.
    E1p=mp*(.5*s*(gamma_e-1.)+\
        np.sqrt(1.+.5*(1.+s*s)*(gamma_e-1.)+\
        .25*s*s*(gamma_e-1.)**(2.))-1.)
    #integrand is the cross section times the proton spectrum
    integrand= lambda Ep: fun.log_interp(vec_p,NEp,Ep)*KO_xsec(Ee,Ep)
    
    #values where the integral is split up
    endpts=[7e-2,9e-2,1e-1,1.5e-1,2e-1,3e-1,4e-1,7e-1,1e0,2e0,4e0,1e1,2e1,\
            4e1,8e1,2e2,8e2] 

    #integration over proton energy
    if E1p<2.5:
        return 1616.5390*density*fun.segmented_int(integrand,endpts,E1p,Epmax)
    else: #large enough E1p doesn't need segmented method
        return 1616.5390*density*fun.fint(integrand,E1p,Epmax)
    
    
#Outputs secondary electron spectrum values at the frequency values 
#given by the array "x"
def secondary_e_spec(x,density,tau_0):
    pts=proton_spec(density,tau_0)
    vec_p=pts[0]  #array of energy values used in electron spectrum
    NEp=pts[1]  #array of electron spectrum values
    
    vec_KO=[]
    vec_x=[]
    for elt in x: 
        gamma_e=elt/5.1099891e-4 + 1.
        E1p=mp*(.5*s*(gamma_e-1.)+\
                np.sqrt(1.+.5*(1.+s*s)*(gamma_e-1.)+\
                .25*s*s*(gamma_e-1.)**(2.))-1.)
        #ensuring the upper limit is larger than the lower limit
        if Epmax>E1p: 
            vec_KO.append(KO_source(elt,density,vec_p,NEp))
            vec_x.append(elt)
        
    return (vec_x,vec_KO)


##Testing the secondary electorn spectrum
###############################################################################
##creating energy array with points equally spaced in log space and
##deleting first and last point
#Emin=1e-4
#Emax=10**(3.8)
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
