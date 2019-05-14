#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt, exp
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
    nmod = len(jomg)

    nind = 0
    for m in xrange(nmod):
        try:
            ndru = len(jomg[m]['jdru'])
        except:
            ndru = 0
        try:
            nsdr = len(jomg[m]['jsdr'])
        except:
            nsdr = 0
        nper = ndru+2*nsdr+npsd
        nind += nper
    delr = np.zeros(nmod,dtype=float)
    expn = np.zeros(nind,dtype=complex)
    etal = np.zeros(nind,dtype=complex)
    etar = np.zeros(nind,dtype=complex)
    etaa = np.zeros(nind,dtype=float)

    pole, resi, rn, tn = PSD (npsd,BoseFermi=1,pade=pade)

    iind = 0

    for m in xrange(nmod):

        try:
            jdru = jomg[m]['jdru']
        except:
            jdru = []

        try:
            jsdr = jomg[m]['jsdr']
        except:
            jsdr = []

        ndru = len(jdru)
        nsdr = len(jsdr)
    
        for idru in xrange(ndru):

            if len(jdru[idru]) != 2:
                raise ValueError('Invalid drude input')

            lamd, gamd = jdru[idru][0], jdru[idru][1]
            expn[iind] = gamd
            etal[iind] = -2.J*lamd*gamd*fBose(-1.J*gamd/temp,pole,resi,rn,tn,exbe)
            etar[iind] = etal[iind].conj()
            etaa[iind] = abs(etal[iind])
            delr[m] += 2.*lamd*gamd/temp*rn

            iind += 1

        for isdr in xrange(nsdr):

            if len(jsdr[isdr]) != 3:
                raise ValueError('Invalid BO input')

            lams, omgs, gams = jsdr[isdr][0], jsdr[isdr][1], jsdr[isdr][2]
            jind = iind+1
            etaBO = 2.*lams*omgs*omgs*gams
            Delta = omgs*omgs-gams*gams/4.0
            delr[m] +=  etaBO*tn/(temp**3)
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
            if Delta > 0:
                etar[iind] = etal[jind].conj()
                etar[jind] = etal[iind].conj()
                etaa[iind] = etaa[jind] = sqrt(abs(etal[iind])*abs(etal[jind]))
            elif Delta < 0:
                etar[iind] = etal[iind].conj()
                etar[jind] = etal[jind].conj()
                etaa[iind] = abs(etal[iind])
                etaa[jind] = abs(etal[jind])
            else:
                raise ValueError("Not prepared for Delta=0")

            try:
                if amod[m] == 'sdr' or amod[m] == 'all':
                    bdip[iind] = fmod[m]
                    bdip[jind] = fmod[m]
            except:
                pass

            iind += 2

        for ipsd in xrange(npsd):
            zomg = -1.J*pole[ipsd]*temp
            jsmd = sum(jwdru(zomg,x) for x in jdru)
            jsms = sum(jwsdr(zomg,x) for x in jsdr)
            jsum = jsmd+jsms
            expn[iind] = pole[ipsd]*temp            
            etal[iind] = -2.J*resi[ipsd]*temp*jsum   
            etar[iind] = etal[iind].conj()
            etaa[iind] = abs(etal[iind])  

            iind += 1
                            
    return etal, etar, etaa, expn, delr
            
if __name__ == '__main__':

    inidic = {
        "q_cl": 'q',
        "exbe": 0,
        "npsd": 2,
        "pade": 2,
        "temp": 1.0,
        "jomg": [{"jdru":[(0.5,0.5)],"jsdr":[(0.1,1.0,0.1)]},{"jdru":[(0.1,0.1)],"jsdr":[(0.1,1.0,0.1)]}],
	"modeFile": "inp_mode.mat",
	"etalFile": "inp_etal.mat",
	"etarFile": "inp_etar.mat",
	"etaaFile": "inp_etaa.mat",
	"expnFile": "inp_expn.mat",
	"delrFile": "inp_delr.mat"
    }
    init(inidic)
    with open('input.json','w') as f:
        json.dump(inidic,f,indent=4) 
    
