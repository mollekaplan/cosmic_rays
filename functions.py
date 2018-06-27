#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:14:31 2018

@author: mollykaplan

Functions needed to find the final proton and electron spectra.
 
"""
import numpy as np
import scipy.integrate as integrate


#the integration function. f is the integrand, Elow is the
#lower limit, Ehi is the upper limit, expb is the exponent on E
#in b(E), and expt is the exponent on E in tau(E)    
def fint(f,Elow,Ehi): 
    #add "points" so quad runs over any range in [10^-4,10^5.7]
    return integrate.quad(f,Elow,Ehi,\
           points=[1e-3,1e-2,1e-1,0.28172800,1e0,1e1,10**(1.5),1e2,10**(2.5),\
           1e3,10**(3.5),1e4,10**(4.5),1e5,10**(5.4)])[0]
    
    
#green function, where E is the lower bound and Ep is the upper bound
def green(f,b,E,Ep):
    int_val= fint(f,E,Ep)
    return 1./(b(E))*np.exp(-int_val)


#injecection spectrum. K is the constant and expo is the exponent on Ep 
def inSpec(K,s,Ep):
    return K*Ep**(s)


#the number spectrum, with E as the lower limit of integration and
#Emax as the upperlimit
def Nspec(f,b,E,Emax,K,s):
    integrand=lambda Ep: inSpec(K,s,Ep)*green(f,b,E,Ep)
    return fint(integrand,E,Emax)


#function to find closest point to some value in an array. returns index.
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#function to find y(val) using log interpolation on y(x), where x and
#y are both arrays.
def log_interp(x,y,val,base10=False):
    i=find_nearest(x,val) #index of array value closest to val
    if x[i]>val and i>0: i-=1  #make sure value is in between x[i] and x[i+1]
    if i==len(x)-1: i-=1  #make sure x[i+1] exists
    
    if not base10:
        slope=(np.log(y[i])-np.log(y[i+1]))/ \
              (np.log(x[i])-np.log(x[i+1]))      
        intercept=np.log(y[i])-slope*np.log(x[i])
    
        #find point on line at val, then exponentiate to return to linear space
        return np.exp(slope*np.log(val)+intercept)
    
    else:
        slope=(np.log10(y[i])-np.log10(y[i+1]))/ \
              (np.log10(x[i])-np.log10(x[i+1]))      
        intercept=np.log10(y[i])-slope*np.log10(x[i])
    
        #find point on line at val, then exponentiate to return to linear space
        return 10**(slope*np.log10(val)+intercept)


#Integration function for more complex integrands, usually ones which
#require interpolation. f is the integrand, endpts is an array of values
#where the integration is broken up, Elow is the lower limit of integration, 
#and Ehi is the upper limit.
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
    