#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 12:08:14 2018

@author: mollykaplan

Parameterizations for the secondary lepton spectrum from charged pion
decay from Kamae 2006. Includes functions to find the diffractive, 
non-diffractive and resonant cross sections for electrons and positrons.
"""
import numpy as np

"""
Helper function to calculate the paramaterization function for some cross 
sections.
"""
def F_array(nel,F_1,F_2,b_0,b_4):
    F=np.zeros(nel)
    for bb in range(nel): 
        if F_1[bb] < -27 and F_2[bb] < -27: F[bb]=0.
        elif F_1[bb] < -27 and F_2[bb] > -27: F[bb]=b_4 * np.exp(F_2[bb])
        elif F_1[bb] > -27 and F_2[bb] < -27: F[bb]=b_0 * np.exp(F_1[bb])
        else: F[bb]=b_0 * np.exp(F_1[bb]) + b_4 * np.exp(F_2[bb])
    return F

"""
Calculation of diffractive cross section using the parameterization function F.
"""
def sigma_d(x,Tp,nel,F,ub):
    W_diff=75.
    L_maxd=np.log10(Tp)
    
    W_diffint=W_diff * (x - L_maxd)
    F_kl=np.zeros(len(W_diffint))
    
    for i in range(len(W_diffint)):
        if W_diffint[i] >= 27: F_kl[i]=0.
        elif W_diffint[i] <= ub: F_kl[i]=1.
        else: F_kl[i]= 1./(np.exp(W_diffint[i]) + 1.)
        
        if F_kl[i] < 1e-12: F_kl[i]=0
        
    sigma=F*F_kl
    
    for bb in range(nel): 
        if sigma[bb] < 1e-30: sigma[bb]=0.
#        
    return sigma

"""
Calculation of non-diffractive cross section using the parameterization 
function F, and the parameters W_ND1, W_NDh, L_min and L_max.
"""
def sigma_nd(x,nel,F,W_ND1,W_NDh,L_min,L_max):
    F_nd_kl = (1./(np.exp(W_ND1 * ( L_min - x )) + 1.)) *\
              (1./(np.exp(W_NDh * ( x - L_max)) + 1.))
    
    sigma=F_nd_kl*F
    for bb in range(nel): 
        if sigma[bb] < 1e-30: sigma[bb]=0.
        
    return sigma

"""
Calculation of resonant cross section using the parameterization 
functions F_kl, d_0, d_1, d_2, d_3 and d_4. The energy bounds 
are lb (lower) and ub (upper).
"""
def sigma_res(x,Tp,nel,F_kl,d_0,d_1,d_2,d_3,d_4):
    F_rexp = -d_1 * ((x-d_2)/(1.0 + d_3*(x-d_2) + d_4*(x-d_2)**2.))**2.    
    F_res2=np.zeros(nel)
    
    for bb in range(nel):
        if F_rexp[bb] > -27: F_res2[bb]=d_0*np.exp(F_rexp[bb])
    
    sigma=F_res2*F_kl
    for bb in range(nel): 
        if sigma[bb] < 1e-30: sigma[bb]=0.
    
    return sigma



############### Functions for calculating cross sections ######################
"""
Calculation of electron diffractive cross section. Valid for energies 
1.94 GeV to 512 TeV. The electron energies are given by the numpy array
x, and the proton energy is given by Tp. nel is the number of points in
the electron spectrum.
"""
def e_d(x,Tp,nel):
    if Tp < 1.95: 
        sigma = np.zeros(nel)
        return sigma
    
    else:
        y = np.log10(Tp) - 3.
    
        ########### diffractive cross section ############
    
        z1 = y + 1.6878
        z2 = y + 9.6400
        b_0=0.20463*np.tanh(-6.2370*(y + 2.2)) - 0.16362*z1*z1 + \
            3.5183e-4*z2**4.
        
        pow1=(y + 2.0154)/(1.0 + 0.62779*(y + 2.0154))
        b_1=1.6537 if -3.2027*pow1*pow1 < -27 \
            else 1.6537 + 3.8530*np.exp(-3.2027*pow1*pow1)
        
        z3 = y + 256.63
        b_2=-10.722 + 0.082672*np.tanh(1.8879*(y + 2.1)) + 0.00014895*z3*z3
        
        pow2 = (y + 1.9877)/(1.0 + 0.40300*(y + 1.988))
        b_3= -0.023752 if -3.3087*pow2*pow2 < -27 \
            else -0.023752 - 0.51734*np.exp(-3.3087*pow2*pow2)
        
        if Tp < 5.51: b_0,b_1,b_2,b_3=0.,0.,0.,0.
        
        z4 = y + 2.9
        b_4=0.94921 + 0.12280*z4*z4 - 7.1585e-4*z4**4. + 0.52130*np.log10(z4)
        b_5=-4.2295 - 1.0025*np.tanh(9.0733*(y + 1.9)) - 0.11452*(y - 62.382)
        b_6=1.4862 + y*(0.99544 + y*(-0.042763 + y*(-0.0040065 + 0.0057987*y)))
        
        z5 = y - 2.8542
        b_7=6.2629 + 6.9517*np.tanh(-0.36480*(y + 2.1)) - 0.26033*z5*z5
        
        F_1=-b_1 * ((x - b_2)/(1. + b_3 * (x - b_2)))**2.
        F_2=-b_5 * ((x - b_6)/(1. + b_7 * (x - b_6)))**2.
        F_diff=F_array(nel,F_1,F_2,b_0,b_4)
        
        ###################    weight function    ###################
        
        return sigma_d(x,Tp,nel,F_diff,-15)


"""
Calculation of electron non-diffractive cross section. The electron energies 
are given by the numpy array x, and the proton energy is given by Tp. 
nel is the number of points in the electron spectrum.
"""
def e_nd(x,Tp,nel):
    y = np.log10(Tp) - 3.

    ########### non-diffractive cross section ############

    z = y + 3.3
    a_0=z*(-0.018639 + z*(2.4315 + z*(-0.57719 + 0.063435*z)))
    a_1=7.1827e-6 + y*(-3.5067e-6 + y*(1.3264e-6 + y*(-3.3481e-7 + \
                                    y*(2.3551e-8 + 3.4297e-9*y))))
    
    z1 = y + 7.9031
    a_2=563.91 - 362.18*np.log10(2.7187*(y + 3.4)) - 2.8924e4/(z1*z1)
    a_3=0.52684 + y*(0.57717 + y*(0.0045336 - 0.0089066*y))
    
    z2 = y + 3.32
    a_4=z2*(0.36108 + z2*(1.6963 + z2*(-0.074456 + z2*(-0.071455 + 0.010473*z2))))
    a_5=9.7387e-5 + 7.8573e-5*np.log10(0.0036055*(y + 4.3)) + \
        0.00024660/(y + 4.9390) - 3.8097e-7*y*y
    a_6=-273.00 - 106.22*np.log10(0.34100*(y + 3.4)) + 89.037*y - 12.546*y*y
    
    z3 = y + 8.5518
    a_7=432.53 - 883.99*np.log10(0.19737*(y + 3.9)) - 4.1938e4/(z3*z3)
    a_8=-0.12756 + y*(0.43478 + y*(-0.0027797 - 0.0083074*y))
    
    rpow=(y + 3.26)/(1. + 9.21*(y + 3.26))
    r_yl=3.63*np.exp(-106*rpow*rpow) + y*(-0.182 - 0.175*y)
    r_yg=1.01
    
    F_tem1=-a_1*((x - a_3) + a_2*(x - a_3)**2.)**2.
    F_tem2=-a_5*((x - a_8) + a_6*(x - a_8)**2. + a_7*(x - a_8)**3.)**2.
    F_nd = a_0*np.exp(F_tem1) + a_4*np.exp(F_tem2)
    
    thresh=15.6
    
    F_nd_fin = r_yl * F_nd if Tp < thresh else r_yg * F_nd
    
    
    ###################    weight function    ###################
    
    W_ND1=20.
    W_NDh=45.
    L_max=0.96*np.log10(Tp)
    L_min=-2.6
    
    return sigma_nd(x,nel,F_nd_fin,W_ND1,W_NDh,L_min,L_max)


"""
Calculation of electron resonant cross section. Valid for energies 0.69 GeV 
to 2.76 GeV. The electron energies are given by the numpy array x, and the 
proton energy is given by Tp. nel is the number of points in the electron 
spectrum.
"""
def e_res(x,Tp,nel):
    if Tp < 0.69 or Tp > 2.76: 
        sigma=np.zeros(nel)
        return sigma 
    
    else:
        y = np.log10(Tp) - 3.
    
        ################### resonance-1600 ####################
    
        pow1=(y + 2.9537)/(1.0 + 1.5221*(y + 2.9537))
        d_0= -0.059458 + 0.0096583*y*y if -56.826*pow1*pow1 < -27 \
             else 0.37790*np.exp(-56.826*pow1*pow1) - 0.059458 + 0.0096583*y*y
        d_1=-5.5135 - 3.3988*y
        d_2=-7.1209 - 7.1850*np.tanh(30.801*(y + 2.1)) + 0.35108*y
        d_3=-6.7841 - 4.8385*y - 0.91523*y*y
        d_4=-134.03 - 139.63*y - 48.316*y*y - 5.5526*y*y*y
        
        
        ##################### energy conservation #####################
        
        W_diff = 75.0
        L_maxr = y + 3.0
        
        w_exp = W_diff*(x - L_maxr)
        
        rexp2=np.where(w_exp <= -27)[0] 
        rexp3=np.intersect1d(np.where(w_exp < 27)[0], np.where(w_exp > -27)[0])    
        F_kl=np.zeros(nel)
        F_kl[rexp2]=1.
        F_kl[rexp3]=1./(np.exp(w_exp[rexp3]) + 1.)
        
        ################ get the resonance cross section ###############
    
        return sigma_res(x,Tp,nel,F_kl,d_0,d_1,d_2,d_3,d_4)

"""
Calculation of positron diffractive cross section. Valid for energies 
1.94 GeV to 512 TeV. The positron energies are given by the numpy array
x, and the proton energy is given by Tp. nel is the number of points in
the positron spectrum.
"""
def p_d(x,Tp,nel):
    if Tp < 1.95: 
        sigma = np.zeros(nel)
        return sigma
    
    else:
        y = np.log10(Tp) - 3.
        
        ########### diffractive cross section ############
        
        z1 = y + 0.67500
        z2 = y + 9.0824
        b_0=(29.192*np.tanh(-0.37879*(y + 2.2)) - 3.2196*z1*z1 + 3.6687e-3*z2**4.)
        pow1 = (y + 1.8781)/(1.0 + 3.8389*(y + 1.8781))
        b_1 = -142.97 if -0.37194*pow1*pow1 < -27 \
              else -142.97 + 147.86*np.exp(-0.37194*pow1*pow1)
        
        z3 = y + 234.65
        b_2=-14.487 - 4.2223*np.tanh(-13.546*(y + 2.2)) + 0.00016988*z3*z3
        pow2 = (y + 1.8194)/(1.0 + 0.99946*(y + 1.8194))
        b_3 = -0.0036974 if -6.1527*pow2*pow2 < -27 \
              else -0.0036974 - 0.41976*np.exp(-6.1527*pow2*pow2)
        
        if Tp < 11.: b_0,b_1,b_2,b_3=0.,0.,0.,0.
        
        z4 = y + 2.95
        pow3 = y + 2.29 - 0.18967*(y + 2.29)
        b4rest=1.8108 + z4*z4*(0.18545 - 2.0049e-3*z4*z4)
        b_4 = b4rest if -14.987*pow3*pow3 < -27 \
              else b4rest + 0.85084*np.exp(-14.987*pow3*pow3)
        b_5=(2.0404 - 0.51548*np.tanh(2.2758*(y + 1.9)) - 0.035009*(y - 6.6555))
        b_6=1.5258 + y*(1.0132 + y*(-0.064388 + y*(-0.0040209 + 0.0082772*y)))
        
        z5 = y - 2.7718
        b_7=3.0551 + 3.5240*np.tanh(-0.36739*(y + 2.1)) - 0.13382*z5*z5
        
        F_1=-b_1 * ((x-b_2)/(1. + b_3 * (x - b_2)))**2.
        F_2=-b_5 * ((x - b_6)/(1. + b_7 * (x - b_6)))**2.
        F_diff=F_array(nel,F_1,F_2,b_0,b_4)
        
        ###################    weight function    ###################
        
        return sigma_d(x,Tp,nel,F_diff,-27)



"""
Calculation of positron non-diffractive cross section. The positron energies 
are given by the numpy array x, and the proton energy is given by Tp. 
nel is the number of points in the positron spectrum.
"""
def p_nd(x,Tp,nel):
    y = np.log10(Tp) - 3.

    ########### non-diffractive cross section ############

    z = (y + 3.3)
    a_0=(z*(-0.79606 + z*(7.7496 + z*(-3.9326 + z*(0.80202 - 0.054994*z)))))
    a_1=(6.7943e-6 + y*(-3.5345e-6 + y*(6.0927e-7 + \
        y*(2.0219e-7 + y*(5.1005e-8 - 4.2622e-8*y)))))
    a_2=(44.827 + 81.378*np.log10(0.027733*(y + 3.5)) - \
         1.3886e4/((y + 8.4417)*(y + 8.4417)))
    a_3=(0.52010 + y*(0.59336 + y*(0.012032 - 0.0064242*y)))
    
    z1 = (y + 3.32)
    a_4=(z1*(2.1361 + z1*(1.8514 + z1*(-0.47872 + z1*(0.0032043 + 0.0082955*z1)))))
    a_5=(1.0845e-6 + 1.4336e-6*np.log10(0.0077255*(y + 4.3)) + \
         0.00013018/((y + 4.8188)*(y + 4.8188)) + 9.3601e-8*y)
    a_6=(-267.74 + 14.175*np.log10(0.35391*(y + 3.4)) + y*(64.669 - 7.7036*y))
    a_7=(138.26 - 529.84*np.log10(0.12467*(y + 3.9)) - \
         1.9869e4/((y + 7.6884)*(y + 7.6884)) + 1.0675*y*y)
    a_8=(-0.14707 + y*(0.40135 + y*(0.0039899 - 0.0016602*y)))
    
    pow1 = (y + 3.25)/(1 + 10.4*(y + 3.25))
    r_yl=0. if -98.9*pow1*pow1 < -27 else 2.22*np.exp(-98.9*pow1*pow1)
    
    F_tem1=-a_1 * ((x - a_3) + a_2 * (x - a_3)**2.)**2.
    F_tem2=-a_5 * ((x - a_8) + a_6 * (x - a_8)**2. + a_7 * (x - a_8)**3.)**2.
    F_nd=F_array(nel,F_tem1,F_tem2,a_0,a_4)    

    thresh=5.52
    F_nd_fin = r_yl * F_nd if Tp < thresh else F_nd
    
    
    ###################    weight function    ###################

    W_ND1=15.
    W_NDh=47.
    L_max=0.94 * np.log10(Tp)
    L_min=-2.6
    
    return sigma_nd(x,nel,F_nd_fin,W_ND1,W_NDh,L_min,L_max)
       

"""
Calculation of positron resonant cross section. The positron energies 
are given by the numpy array x, and the proton energy is given by Tp. 
nel is the number of points in the positron spectrum.

Reso 1 is valid for energies 0.488 GeV to 1.95 GeV. Reso 2 is valid
for energies 0.69 GeV to 2.76  GeV.
"""
def p_res(x,Tp,nel,reso=1):
    if reso==1: lb,ub=0.488,1.95
    else: lb,ub=0.69,2.76
    
    if Tp < lb or Tp > ub: 
        sigma=np.zeros(nel)
        return sigma 
    else:
        y = np.log10(Tp) - 3.
    
        ################### resonance-1600 ####################
    
        if reso==1:
            if y==0: y=0.1
    
            pow1 = ((y + 3.1272)/(1.0 + 0.22831*(y + 3.1272)))
            rest= -(6.5855 + 9.6984/y - 0.41256*y*y)
            d_0=rest if -67.857*pow1*pow1 < -27 else 2.9841*np.exp(-67.857*pow1*pow1) + rest
            d_1=6.8276 + 5.2236*y + 1.4630*y*y
            d_2=-6.0291 - 6.4581*np.tanh(5.0830*(y + 2.1)) + 0.46352*y
            d_3=0.59300 + 0.36093*y
            d_4=0.77368 + 0.44776*y + 0.056409*y*y
        
        else:
            pow1 = ((y + 2.9485)/(1.0 + 1.2892*(y + 2.9485)))
            rest=-(0.23720 - 0.041315*y*y)
            d_0 = rest if -56.544*pow1*pow1 < -27 else 1.9186*np.exp(-56.544*pow1*pow1) + rest
            d_1=-4.9866 - 3.1435*y
            d_2=-7.0550 - 7.2165*np.tanh(31.033*(y + 2.1)) + 0.38541*y
            d_3=-2.8915 - 2.1495*y - 0.45006*y*y
            d_4=-1.2970 - 0.13947*y + 0.41197*y*y + 0.10641*y*y*y
    
        
        ##################### energy conservation #####################
        
        W_diff = 75.0
        L_maxr = y + 3.0
        
        w_exp = W_diff*(x - L_maxr)
        
        F_kl=1./(np.exp(w_exp) + 1.)
    
        ################ get the resonance cross section ###############
        
        if reso==1:
            return sigma_res(x,Tp,nel,F_kl,d_0,d_1,d_2,d_3,d_4)
        else:
            return sigma_res(x,Tp,nel,F_kl,d_0,d_1,d_2,d_3,d_4)
    
    
    


    
