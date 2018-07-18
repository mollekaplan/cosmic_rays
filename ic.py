#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 10:30:14 2018

@author: mollykaplan

Calulates the electron energy loss rate from inverse compton losses,
using an input IR spectrum (here, uses a model for NGC 253, taken
from ir.csv). Loads in the integral over q, which is calulated off
line in qint.py. 

b_ic returns an array of loss values at a given array of electron
energy values. These can then be read into the electrons.py, so
their contribution to the energy loss rate can be considered.

Under testing section, can output a plot of the b_ic at an array
of energy electron values.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate
import functions as fun

#importing spectral data
data_file='ir.csv'
col_names=['wl','flux']
data = pd.read_csv(data_file,names=col_names)
qint=np.load("qint.npy")
qint_E=np.load("qint_E.npy")
mUph=np.load("mUph.npy")

"""
Function to create arrays out of spectral data, returns tuple with
IR spectrum energy and spectral data.
"""
def spec_data(Uph):
    wl=data.wl
    spec=data.flux
    
    #turning wavelengths into energy values in eV
    e_ir=[]
    for elt in wl: e_ir.append(1.2431/elt)
    
    #converting the spectrum into required units
    ne=[]
    for i in range(len(spec)):
        ne.append(125.85*Uph/mUph*spec[i]/e_ir[i])
    
    return (e_ir,ne)
   
""" 
Function which outputs the IC energy loss rate at an array of 
energy value given by E. Thompson limit energy from Schlickeiser
p.89. The full integral is taken from Schlickeiser p.88.
"""      
def b_ic(E,Uph):
    #read spectrum data from spec_data
    e_ir=spec_data(Uph)[0]
    ne=spec_data(Uph)[1]
    
    #log of IR energies
    l_e_ir=[]
    for elt in e_ir: l_e_ir.append(np.log10(elt))
    
    ind=0 #index to find changeover from Thompson to KN
    b_arr=[] #empty array for storing b_ic values
    co=len(E)-1 #start at highest index

    for j in range(len(E)): 
        #Thompson limit:
        if E[j]<1e-1: 
            gam=E[j]/5.1099891e-4 + 1.
            b_arr.append(8.36e-10*Uph*gam*gam)
            
        #KN limit    
        else: 
            if ind<co: co=ind #find changeover index
            
            #find electron energy that labels the qint array
            k=fun.find_nearest(qint_E,E[j])
            q_arr=qint[k]
            
            #fill array with integrand values
            int_arr=[]
            for i in range(len(e_ir)): 
                int_arr.append(ne[i]/e_ir[i]*q_arr[i])           

            #integrate int_arr over IR spectrum energies, ignore
            #energies above the 279th index since the spectrum is
            #not well behaved above these values
            b_arr.append(123.2*integrate.simps(int_arr[0:280],e_ir[0:280]))
        
        ind+=1
    
    #apply conversion to match up Thompson and KN results
    conv=8.36e-10*Uph*(E[co]/5.1099891e-4 + 1.)**2./b_arr[co]
    for i in range(co,len(b_arr)):
        b_arr[i]*=conv
    
    return b_arr


###Testing the b_ic values
###############################################################################
#Uph=1e4
#Emin=1e-4
#Emax=1e6
#npts=100
#E=np.geomspace(Emin,Emax,npts)
#
#e_ir=spec_data(Uph)[0]
#
#arr=b_ic(E,Uph)
#
#print(arr)
#
##plotting the number spectrum
#f1 = plt.figure(figsize=(11,6))
#ax1 = f1.add_subplot(111)
#plt.plot(E,arr)
#ax1.set_xscale('log')
#ax1.set_yscale('log')

