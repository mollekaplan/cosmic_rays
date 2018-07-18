#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 11:24:22 2018

@author: mollykaplan

Calculates the q integral in the inverse Compton losses cross 
section, as given Schlickeiser p.88. Saves the results (along with 
the energy values chosen) so this only needs to be run once.
"""
import numpy as np
import pandas as pd
import functions as fun

#read spectrum data and create arrays
data_file='ir.csv'
col_names=['wl','flux']
data = pd.read_csv(data_file,names=col_names)
wl=data.wl
e_ir=[]
for elt in wl: e_ir.append(1.2431/elt)

#electron energy values
Emin=1e-4
Emax=1e6
npts=100
E=np.geomspace(Emin,Emax,npts)

"""
Outputs array of values for the integral over q, each with a 
different ir spectrum energy. 
"""
def qint(E,e):
    gam=E/5.1099891e-4 + 1.    
    q_arr=[]
    for elt in e: 
        ge=4.*elt*gam/5.1099891e5
        integrand= lambda q: ge*ge*q/(1.+ge*q)**3.*(2.*q*np.log(q) + \
                        (1.+2.*q)*(1.-q) + \
                        0.5*(ge*q)**2.*(1.-q)/(1.+ge*q))
        q_arr.append(fun.fint(integrand,0.,1.))
    
    return q_arr

"""
Create an array of qint values, where each point in the matrix
represents a different electron energy and ir spectrum energy
pair
"""
full_int=np.array(qint(E[0],e_ir))
for i in range(1,len(E)):
    full_int=np.vstack([full_int,qint(E[i],e_ir)])

#save the results
np.save("qint_E",E)
np.save("qint",full_int)