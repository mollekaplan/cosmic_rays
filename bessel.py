#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 10:42:08 2018

@author: mollykaplan

Calculates F(x)=x*int_x^(infinity) K_(5/3)(y) dy for various values 
of x. The function Fx outputs a tuple with the x values in the first
index and F(x) values in the second.
"""

import numpy as np
from scipy.special import kv
import scipy.integrate as integrate

#creating F(x) function
def mod_bessel(x):
    return x*integrate.quad(lambda xi: kv(5/3,xi), x, float("inf"))[0]


#range of values for x, equally spaced in log space
x=np.geomspace(1e-5,1e2,100)

#creating an array of values from integration of the 
#modified bessel function
F=[]
for elt in x: F.append(mod_bessel(elt))

#outputting a tuple with arrays of (1) x values and (2) F values
def Fx():
    return (x,F)
