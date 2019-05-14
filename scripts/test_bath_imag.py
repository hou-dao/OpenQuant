#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt, exp, pi
from BoseFermiExpansion import PSD
import armadillo as arma 
from cmath import exp as cexp
import json

def jwdru (omg,jdru):
    lamd, gamd = jdru[0], jdru[1]
    return 2.0*lamd*gamd*omg/(omg**2+gamd**2)

def jwsdr (omg,jsdr):
    lams, omgs, gams = jsdr[0], jsdr[1], jsdr[2]
    return 2.0*lams*omgs**2*gams*omg/((omg**2-omgs**2)**2+(omg*gams)**2)

def fBose (x,pole,resi,rn,tn,sign=0):
    if sign == 0:
        return 1/x+0.5+rn*x+tn*x**3+sum(2.0*resi[i]*x/(x**2+pole[i]**2) for i in xrange(len(pole)))
    else:
        if isinstance(x,complex):
            return 1.0/(1-cexp(-x))
        else:
            return 1.0/(1-exp(-x))
    
def generate (temp,npsd,pade,jomg):
        
    exbe = 0 
    nind = 0
    ndru = 0
    nsdr = 0
    if (jomg.has_key('jdru')):
        ndru = 1
    if (jomg.has_key('jsdr')):
        nsdr = 1
    nind = ndru+2*nsdr+npsd
    expn = np.zeros(nind,dtype=complex)
    etal = np.zeros(nind,dtype=complex)

    pole, resi, rn, tn = PSD (npsd,BoseFermi=1,pade=pade)

    iind = 0
    
    if ndru == 1:
        jdru = jomg['jdru']
        if len(jdru) != 2:
            raise ValueError('Invalid drude input')
        lamd, gamd = jdru[0], jdru[1]
        expn[iind] = gamd
        etal[iind] = -2.J*lamd*gamd*fBose(-1.J*gamd/temp,pole,resi,rn,tn,exbe)
        iind += 1

    if nsdr == 1:
        jsdr = jomg['jsdr']
        if len(jsdr) != 3:
            raise ValueError('Invalid BO input')
        lams, omgs, gams = jsdr[0], jsdr[1], jsdr[2]
        jind = iind+1
        etaBO = 2.*lams*omgs*omgs*gams
        Delta = omgs*omgs-gams*gams/4.0
        if Delta > 0:
            OmgB = sqrt(Delta)
            expn[iind] = 0.5*gams+1.J*OmgB
            expn[jind] = 0.5*gams-1.J*OmgB
        elif Delta < 0:
            OmgB = sqrt(-Delta)
            expn[iind] = 0.5*gams+OmgB
            expn[jind] = 0.5*gams-OmgB           
        else:
            raise ValueError("Not prepared for Delta=0")
        z1, z2 = -1.J*expn[iind], -1.J*expn[jind]
        etal[iind] = -2.J*etaBO*z1/(2.*z1*(z1+z2)*(z1-z2))*fBose(z1/temp,pole,resi,rn,tn,exbe)
        etal[jind] = -2.J*etaBO*z2/(2.*z2*(z2+z1)*(z2-z1))*fBose(z2/temp,pole,resi,rn,tn,exbe)
        iind += 2

    for ipsd in xrange(npsd):
        zomg = -1.J*pole[ipsd]*temp
        jsum = 0.0+0.0J
        if ndru == 1:
            jsum += jwdru(zomg,jomg['jdru'])
        if nsdr == 1:
            jsum += jwsdr(zomg,jomg['jsdr'])
        expn[iind] = pole[ipsd]*temp            
        etal[iind] = -2.J*resi[ipsd]*temp*jsum   
        iind += 1
                        
    return etal, expn
            
if __name__ == '__main__':

    temp = 1.0
    pade = 2
    jomg = {'jdru':(1.0,5.0)}

    etas1, gams1 = generate(temp,2,pade,jomg)
    etas2, gams2 = generate(temp,16,pade,jomg)
    etas3, gams3 = generate(temp,64,pade,jomg)
    etas4, gams4 = generate(temp,256,pade,jomg)
    etas5, gams5 = generate(temp,1024,pade,jomg)

    with open('ctau.dat','w') as f:
        dt =1.0/(temp*1000)
        for i in xrange(1000):
            ct1 = np.dot(etas1,np.exp(1.J*gams1*i*dt))
            ct2 = np.dot(etas2,np.exp(1.J*gams2*i*dt))
            ct3 = np.dot(etas3,np.exp(1.J*gams3*i*dt))
            ct4 = np.dot(etas4,np.exp(1.J*gams4*i*dt))
            ct5 = np.dot(etas5,np.exp(1.J*gams5*i*dt))
            print >> f, 11*'%16.6e'%(i*dt, \
                                     ct1.real,ct1.imag, \
                                     ct2.real,ct2.imag, \
                                     ct3.real,ct3.imag, \
                                     ct4.real,ct4.imag, \
                                     ct5.real,ct5.imag)
