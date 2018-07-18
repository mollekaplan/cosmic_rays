#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 10:45:45 2018

@author: mollykaplan

Calculates the total flux over the spectrum in ir.csv (which uses
data for NGC 253). Saves total flux.
"""

import numpy as np
import functions as fun
import pandas as pd
from scipy.interpolate import interp1d

#load data
data_file='ir.csv'
col_names=['wl','flux']
data = pd.read_csv(data_file,names=col_names)
wl=data.wl
spec=data.flux

#creating arrays of data
nu_ir=[]
for elt in wl: nu_ir.append(3e14/elt)
l_wl=[]
for elt in nu_ir: l_wl.append(np.log10(elt))
l_spec=[]
for elt in spec: l_spec.append(np.log10(elt))

#only until index 279 to avoid weird regions in the spectrum
integrand=lambda e: 10.**interp1d(l_wl[0:280],l_spec[0:280], 
                                  fill_value=-np.Inf, 
                                  bounds_error=False, 
                                  kind='cubic')(np.log10(e))

#calculate the integral
endpts=[nu_ir[50],nu_ir[100],nu_ir[130],nu_ir[160],nu_ir[190],
        nu_ir[220],nu_ir[250]]
mUph=5.205e-13*(fun.segmented_int(integrand,endpts,nu_ir[0],nu_ir[279]))

#save total flux
np.save("mUph",mUph)