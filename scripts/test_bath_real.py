#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt, exp, pi
from BoseFermiExpansion import PSD
import armadillo as arma 
from cmath import exp as cexp
import json

if __name__ == '__main__':

    # bath
    lamd = 0.2
    gamd = 1.0

    nmod = 1
    ini['bath']['temp'] = temp
    ini['bath']['nmod'] = nmod
    ini['bath']['jomg'] = [{"jdru":[(lamd,gamd)]} for i in xrange(nmod)]
   #ini['bath']['jomg'] = [{"jsdr":[(lams,omgs,gams)]} for i in xrange(nmod)]
    ini['bath']['pade'] = 0
    ini['bath']['npsd'] = 16
    bath.init (ini['bath'])


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
