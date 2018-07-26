#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:14:31 2018

@author: mollykaplan

Functions needed to find the final proton and lepton spectra.
"""
import numpy as np
import scipy.integrate as integrate


"""
An integration function used throughout the analysis. f is the 
integrand, Elow is the lower limit, and Ehi is the upper limit. 
"""
def fint(f,Elow,Ehi): 
    #add "points" so quad runs over any range in [10^-4,10^5.7]
    return integrate.quad(f,Elow,Ehi,
           points=[1e-3,1e-2,1e-1,0.28172800,1e0,1e1,10**(1.5),1e2,10**(2.5),
           1e3,10**(3.5),1e4,10**(4.5),1e5,10**(5.4)])[0]
    
"""
Green function as given by Domingo-Santamaria et. al. 2008, Eq. 2.
E is the lower bound and Ep (E' in the paper) is the upper bound.
f is the integrand, to be input later.
"""
def green(f,b,E,Ep):
    int_val= fint(f,E,Ep)
    return 1./(b(E))*np.exp(-int_val)

"""
Power law injection spectrum. K is the constant and Ep is E', 
to match E' in the Green function. 
"""
def inSpec_p(K,s,Ep):
    return K*(Ep + .9383)**(-s)

"""
The number spectrum as given by Domingo-Santamaria et. al. 2008, 
Eq. 3. E is the lower limit of integration and Emax is the upper 
limit.
"""
def Nspec_p(f,b,E,Emax,K,s):
    integrand=lambda Ep: inSpec_p(K,s,Ep)*green(f,b,E,Ep)
    return fint(integrand,E,Emax)

"""
Injection spectrum for electrons, which includes the secondary 
electron spectrum for knock-ons. Will include secondary spectrum 
from pion decay once that is ready. K is the constant and Ep is E', 
to match E' in the Green function. 
"""
def inSpec_e(K,s,Ep,sec_spec):
    #calculate ratio of leptons to protons, given by Paglione et. al.
    #2012 Eq. 2
    N_ratio=(5.1099891e-4/.9383)**(.5*(s-1.))
    return K*(N_ratio*(Ep + 5.1099891e-4)**(-s) + \
              10.**sec_spec(np.log10(Ep + 5.1099891e-4)))

"""
The number spectrum for electrons, with E as the lower limit of 
integration and Emax as the upper limit.
"""
def Nspec_e(f,b,E,Emax,K,s,sec_spec):
    integrand=lambda Ep: inSpec_e(K,s,Ep,sec_spec)*green(f,b,E,Ep)
    return fint(integrand,E,Emax)

"""
Function to find closest point to some value in an array. 
Returns index.
"""
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

"""
Integration function for more complex integrands, usually ones which
require interpolation. f is the integrand, endpts is an array of 
values where the integration is broken up, Elow is the lower limit 
of integration, and Ehi is the upper limit.
"""
def segmented_int(f,endpts,Elow,Ehi):
    arr=[] #array which will hold integral values
    
    i=find_nearest(endpts,Elow) #finding index of endpts where Elow is closest
    if Elow>endpts[i]: i+=1 #making sure endpt is above min
    if i>=len(endpts): return fint(f,Elow,Ehi) #if no endpts above, return integral
    
    #creating the array
    for j in range(i,len(endpts)-1): 
        arr.append(integrate.quad(f,endpts[j],endpts[j+1])[0])
    
    int_sum=fint(f,Elow,endpts[i]) #first segment
    for j in range(len(arr)): int_sum+=arr[j]  #adding up all other integrals
    int_sum+=fint(f,endpts[-1],Ehi) #final segment
        
    return int_sum
    
