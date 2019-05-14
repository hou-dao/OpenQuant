import sys
sys.path.append('../scripts')

import argparse
import json
import numpy as np
import armadillo as arma
import syst
import bath_generator as bath
from math import sqrt, pi

kt2au = 1.0
cm2au = 1.0
fs2au = 1.0

if __name__ == '__main__':

    ini = {'bath':{},'syst':{},'hidx':{}}

    # bath
    temperature = 1.0
    lamd = 0.5*3.14*5.0
    gamd = 5.0

    nmod = 1
    npsd = 2
    pade = 1
    jomg = [{"jdru":[(lamd,gamd)]} for i in xrange(nmod)]

    etal, etar, etaa, expn, delr = bath.generate(temperature,npsd,pade,jomg)

    temp = np.array([temperature])
    mode = np.array([0 for i in xrange(len(etal))])

    ini['bath']['temp'] = arma.arma2json(temp)
    ini['bath']['mode'] = arma.arma2json(mode)
    ini['bath']['etal'] = arma.arma2json(etal)
    ini['bath']['etar'] = arma.arma2json(etar)
    ini['bath']['etaa'] = arma.arma2json(etaa)
    ini['bath']['expn'] = arma.arma2json(expn)
    ini['bath']['delr'] = arma.arma2json(delr)

    # syst
    hams = np.zeros((2,2),dtype=complex)
    hams[0,1] = hams[1,0] = 0.5

    qmds = np.zeros((1,2,2),dtype=complex)
    qmds[0,0,0] = 0.5
    qmds[0,1,1] =-0.5

    ini['syst']['ham1'] = arma.arma2json(hams)
    ini['syst']['qmd1'] = arma.arma2json(qmds)
    
    # hidx
    ini['hidx']['trun'] = 0
    ini['hidx']['lmax'] = 20
    ini['hidx']['nmax'] = 10000000
    ini['hidx']['ferr'] = 1.0e-6
    ini['hidx']['expn'] = arma.arma2json(expn)

    sdip = np.zeros((2,2),dtype=complex)
    sdip[0,1] = sdip[1,0] = 1.0

    #proprho
    jsonInit = {"deom":ini,
                "dt": 0.002,
                "nt": 32768,
                "nk": 10,
                "staticErr": 1.0e-6,
                "inistate": 0,
                "sdip1": arma.arma2json(sdip),
                "sdip2": arma.arma2json(sdip)
            }

    json.encoder.FLOAT_REPR = lambda f: ("%16.6e"%f)
    with open('input.json','w') as f:
        json.dump(jsonInit,f,indent=4) 
