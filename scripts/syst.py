#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from unitconverter import cm2au
import armadillo as arma

def init (inidic,hams,qmds):
        
    arma.save (hams,inidic['hamsFile'])
    arma.save (qmds,inidic['qmdsFile'])
