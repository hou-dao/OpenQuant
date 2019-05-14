import json
import numpy as np
import armadillo as arma
import syst
import bath

if __name__ == '__main__':

    with open('_default.json') as f:
        ini = json.load(f)
    
    # syst
    hams = np.zeros((2,2),dtype=complex)
    hams[0,0] = hams[0,1] = hams[1,0] = 1.0
    hams[1,1] = -1.0
    qmds = np.zeros((1,2,2),dtype=complex)
    qmds[0,0,0] =  1.0
    qmds[0,1,1] = -1.0
    dips = np.zeros((2,2),dtype=complex)
    syst.init (ini['syst'],hams,qmds,dips)
   
    # bath
    ini['bath']['temp'] = 1.0
    ini['bath']['jsdr'] = [{'lams':1.0,'omgs':5.0,'gams':1.0}]
    ini['bath']['nmod'] = qmds.shape[0]
    ini['bath']['pade'] = 2
    ini['bath']['npsd'] = 1
    bath.init (ini['bath'])
    
    # hidx
    ini['hidx']['lmax'] = 15
    ini['hidx']['nmax'] = 30000
    ini['hidx']['ferr'] = 2.0e-6
    
    with open('input.json','w') as f:
        json.dump(ini,f,indent=4) 
