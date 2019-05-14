#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt, exp
from BoseFermiExpansion import PSD
import armadillo as arma 
from cmath import exp as cexp
import json

# Drude
def jwdru (omg,jdru):
    lamd, gamd = jdru[0], jdru[1]
    return 2.0*lamd*gamd*omg/(omg**2+gamd**2)

# Super-Drude
def jwsdr (omg,jsdr):
    lams, omgs, gams = jsdr[0], jsdr[1], jsdr[2]
    return 2.0*lams*omgs**2*gams*omg/((omg**2-omgs**2)**2+(omg*gams)**2)

# Quadratic case
def jwqdr (omg,jqdr):
    lamq, omgq, lamp, gamp = jqdr[0], jqdr[1], jqdr[2], jqdr[3]
    tmp1 = 2*lamq*omgq
    tmp2 = lamp*gamp
    tmp3 = omg**2-omgq**2
    return tmp1*tmp2*omgq*gamp*omg/(((tmp3-tmp2)*omg)**2+(gamp*tmp3)**2)

def qdrgm (jqdr):
    lamq, omgq, lamp, gamp = jqdr[0], jqdr[1], jqdr[2], jqdr[3]
    qn = np.zeros(7,dtype=float)
    qn[0] = 1.0
    qn[2] = gamp**2-2.0*(omgq**2+lamp*gamp)
    qn[4] = (omgq**2+lamp*gamp)**2-2*(omgq*gamp)**2
    qn[6] = (omgq*omgq*gamp)**2
    rts = np.roots(qn)
    gms = 1.J*rts[np.where(rts.imag<0)]
   #print 'rts:', rts
   #print 'gms:', gms
    i0, i1, i2 = -1, -1, -1
    if (all(abs(gms.imag)<1.e-10)): # All pure imaginary pole
        i0, i1, i2 = 0, 1, 2
        gms = gms.real
        return gms[0], gms[1], gms[2]
    elif abs(gms[0].imag)<1.e-10 and abs(gms[1].imag+gms[2].imag)<1.e-10:
        i0, i1, i2 = 0, 1, 2
    elif abs(gms[2].imag)<1.e-10 and abs(gms[0].imag+gms[1].imag)<1.e-10:
        i0, i1, i2 = 2, 0, 1
    else:
        i0, i1, i2 = 1, 0, 2
    if abs(gms[i1].real-gms[i2].real)>1.e-10:
        raise ValueError('Hey buddy, you got a problem!')
    re = 0.5*(gms[i1].real+gms[i2].real)
    im = 0.5*(gms[i1].imag-gms[i2].imag)
    gm0 = gms[i0].real
    gm1 = re+1.J*im
    gm2 = re-1.J*im
    return gm0, gm1, gm2

def fBose (x,pole,resi,rn,tn,sign=0):
    if sign == 0:
        return 1/x+0.5+rn*x+tn*x**3+sum(2.0*resi[i]*x/(x**2+pole[i]**2) for i in xrange(len(pole)))
    else:
        if isinstance(x,complex):
            return 1.0/(1-cexp(-x))
        else:
            return 1.0/(1-exp(-x))
    
def init (inidic):
        
    try:
        q_cl = inidic['q_cl']
        exbe = inidic['exbe']
    except:
        q_cl = 'q'
        exbe = 0

    npsd = inidic['npsd']
    pade = inidic['pade']
    temp = inidic['temp']
    jomg = inidic['jomg']
    nmod = len(jomg)

    nind = 0
    mode = []
    for m in xrange(nmod):
        try:
            ndru = len(jomg[m]['jdru'])
        except:
            ndru = 0
        try:
            nsdr = len(jomg[m]['jsdr'])
        except:
            nsdr = 0
        try:
            nqdr = len(jomg[m]['jqdr'])
        except:
            nqdr = 0
        nper = ndru+2*nsdr+3*nqdr+npsd
        nind += nper
        mode.extend([m for i in xrange(nper)])
    mode = np.array(mode)
    delr = np.zeros(nmod,dtype=float)
    expn = np.zeros(nind,dtype=complex)
    etal = np.zeros(nind,dtype=complex)
    etar = np.zeros(nind,dtype=complex)
    etaa = np.zeros(nind,dtype=float)
    bdip = np.zeros(nind,dtype=float)

    pole, resi, rn, tn = PSD (npsd,BoseFermi=1,pade=pade)

    iind = 0

    for m in xrange(nmod):
        # Dru
        try:
            jdru = jomg[m]['jdru']
        except:
            jdru = []
        # Sdr
        try:
            jsdr = jomg[m]['jsdr']
        except:
            jsdr = []
        # Qdr
        try:
            jqdr = jomg[m]['jqdr']
        except:
            jqdr = []
        # Num
        ndru = len(jdru)
        nsdr = len(jsdr)
        nqdr = len(jqdr)
    
        for idru in xrange(ndru):

            if len(jdru[idru]) != 2:
                raise ValueError('Invalid drude input')

            lamd, gamd = jdru[idru][0], jdru[idru][1]
            expn[iind] = gamd
            etal[iind] = -2.J*lamd*gamd*fBose(-1.J*gamd/temp,pole,resi,rn,tn,exbe)
            if q_cl == 'cl':
                etal[iind] = etal[iind].real
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
                if q_cl == 'cl':
                    etal[iind], etal[jind] = 0.5*(etal[iind]+etal[jind].conj()), \
                            0.5*(etal[jind]+etal[iind].conj())
                etar[iind] = etal[jind].conj()
                etar[jind] = etal[iind].conj()
                etaa[iind] = etaa[jind] = sqrt(abs(etal[iind])*abs(etal[jind]))
            elif Delta < 0:
                if q_cl == 'cl':
                    etal[iind], etal[jind] = etal[iind].real, etal[jind].real
                etar[iind] = etal[iind].conj()
                etar[jind] = etal[jind].conj()
                etaa[iind] = abs(etal[iind])
                etaa[jind] = abs(etal[jind])
            else:
                raise ValueError("Not prepared for Delta=0")

            iind += 2

        for iqdr in xrange(nqdr):

            if len(jqdr[iqdr]) != 4:
                raise ValueError('Invalid Jw input')

            lams, omgs, lamt, gamt = jqdr[iqdr][0], jqdr[iqdr][1], jqdr[iqdr][2], jqdr[iqdr][3]
            jind = iind+1
            kind = iind+2
            expn[iind], expn[jind], expn[kind] = qdrgm(jqdr[iqdr])
            etaBO = 2.0*lams*lamt*(omgs*gamt)**2
            z0, z1, z2 = -1.J*expn[iind], -1.J*expn[jind], -1.J*expn[kind]

            etal[iind] = -1.J*etaBO/((z0**2-z1**2)*(z0**2-z2**2))*fBose(z0/temp,pole,resi,rn,tn,exbe)
            etal[jind] = -1.J*etaBO/((z1**2-z0**2)*(z1**2-z2**2))*fBose(z1/temp,pole,resi,rn,tn,exbe)
            etal[kind] = -1.J*etaBO/((z2**2-z0**2)*(z2**2-z1**2))*fBose(z2/temp,pole,resi,rn,tn,exbe)

            if (all(abs(expn[iind:kind+1].imag)<1.0e-10)):
                etar[iind:kind+1] = etal[iind:kind+1].conj()
                etaa[iind:kind+1] = abs(etal[iind:kind+1])

                if q_cl == 'cl':
                    etal[iind:kind+1] = etal[iind:kind+1].real
                    etar[iind:kind+1] = etal[iind:kind+1].conj()
                    etaa[iind:kind+1] = abs(etal[iind:kind+1])
            else:
                etar[iind] = etal[iind].conj()
                etar[jind] = etal[kind].conj()
                etar[kind] = etal[jind].conj()
                etaa[iind] = abs(etal[iind])
                etaa[jind] = etaa[kind] = sqrt(abs(etal[jind])*abs(etal[kind]))

                if q_cl == 'cl':
                    etal[iind] = etal[iind].real
                    etar[iind] = etal[iind].conj()
                    etaa[iind] = abs(etal[iind])
                    etal[jind], etal[kind] = 0.5*(etal[jind]+etal[kind].conj()), \
                                0.5*(etal[kind]+etal[jind].conj())
                    etar[jind] = etal[kind].conj()
                    etar[kind] = etal[jind].conj()
                    etaa[jind] = etaa[kind] = sqrt(abs(etal[jind])*abs(etal[kind]))

            iind += 3

        for ipsd in xrange(npsd):
            zomg = -1.J*pole[ipsd]*temp
            jsmd = sum(jwdru(zomg,x) for x in jdru)
            jsms = sum(jwsdr(zomg,x) for x in jsdr)
            jsmq = sum(jwqdr(zomg,x) for x in jqdr)
            jsum = jsmd+jsms+jsmq
            expn[iind] = pole[ipsd]*temp            
            etal[iind] = -2.J*resi[ipsd]*temp*jsum   
            etar[iind] = etal[iind].conj()
            etaa[iind] = abs(etal[iind])  

            iind += 1
                            
    arma.save (mode,inidic['modeFile'])
    arma.save (etal,inidic['etalFile'])
    arma.save (etar,inidic['etarFile'])
    arma.save (etaa,inidic['etaaFile'])
    arma.save (expn,inidic['expnFile'])
    arma.save (delr,inidic['delrFile'])
    try:
        arma.save (bdip,inidic['bdipFile'])
    except:
        print 'bdipFile not provided!'
        pass

    return etal, expn
            
            
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
    
