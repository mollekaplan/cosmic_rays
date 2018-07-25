#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 15:32:27 2018

@author: mollykaplan

Uses the parameterization functions in pion_param to create the total
injection spectrum for secondary leptons from charged pion decay.
"""

import pion_param as pp
import numpy as np

"""
Outputs the secondary lepton injection spectrum Q. Ee is an array of 
lepton energy values, nel is the size of Ee, Ep is an array of
proton energies, NEp is the steady state proton spectrum at the 
energies in Ep, and density is the ambient gas density.
"""
def pion_sec(Ee,nel,Ep,NEp,density):
    x=np.log10(Ee)
    
    width = 10.**((4.5 + 2.)/1121.) - 1.
    
    rese=np.zeros(nel)
    resp=np.zeros(nel)
    for bb in range(len(Ep)):
        if Ep[bb]<0.488: continue  #skip low energy protons
        
        #finding electron injections spectrum
        diff_e=pp.e_d(x,Ep[bb],nel)
        Ndiff_e=pp.e_nd(x,Ep[bb],nel)
        r1_e=pp.e_res(x,Ep[bb],nel)
        rese += (NEp[bb]*Ep[bb])*\
                (diff_e + Ndiff_e + r1_e) *width
        
        #finding positron injection spectrum
        diff_p=pp.p_d(x,Ep[bb],nel)
        Ndiff_p=pp.p_nd(x,Ep[bb],nel)
        r1_p=pp.p_res(x,Ep[bb],nel,1)
        r2_p=pp.p_res(x,Ep[bb],nel,2)
        resp += (NEp[bb]*Ep[bb])*\
                (diff_p + Ndiff_p + r1_p + r2_p) * width 
    
    #calultion of total (electron+positron) injection spectrum
    res_tot=rese+resp
    Q=9.4605284e-4*density*res_tot/Ee
    return (Ee,Q)
