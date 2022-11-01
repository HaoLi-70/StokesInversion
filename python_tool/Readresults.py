import sys, os, struct

import numpy as np

################################################################################

def readresult(filename): 
    
    try:
        f = open(filename,'rb')
    except:
        f.close()
        print("Could not open file " + filename)
        return {'check': False}

    try:
        result = {}

        tmp = f.read(4)
        result['Nsamples'] = int(struct.unpack('i', tmp)[0])
        
        tmp = f.read(4)
        result['Num_Par'] = int(struct.unpack('i', tmp)[0])
       
        tmp = f.read(8*result['Nsamples']*result['Num_Par'])
        result['Samples'] = \
            np.reshape(np.array(struct.unpack\
            ('d'*result['Nsamples']*result['Num_Par'], tmp)), \
                       (-1,result['Nsamples']))
        return result
    except:
        print("Could not read file" + filename)
        return {'check': False}

################################################################################
