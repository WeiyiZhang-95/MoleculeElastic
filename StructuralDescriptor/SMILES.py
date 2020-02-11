

import openbabel
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("xyz", "smi")
import numpy as np

import os, sys

def smile(filename, dir=os.getcwd()):
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol,os.path.join(dir,filename)
    outMDL = obConversion.WriteString(mol).split()[0]
    return outMDL
  
  
SMILES_CHARS = [' ',
                  '#', '%', '(', ')', '+', '-', '.', '/',
                  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                  '=', '@',
                  'A', 'B', 'C', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P',
                  'R', 'S', 'T', 'V', 'X', 'Z',
                  '[', '\\', ']',
                  'a', 'b', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'p', 'r', 's',
                  't', 'u']

__smile_max_len__ = 155

smi2index = dict( (c,i) for i,c in enumerate( SMILES_CHARS ) )
def smiles_encoder( smiles, maxlen=__smile_max_len__ ):
    X = np.zeros( ( maxlen, len( SMILES_CHARS ) ) )
    for i, c in enumerate( smiles ):
        if i>=maxlen:
            break
        X[i, smi2index[c] ] = 1
    return X


def test():
    filename = 'linker0-un-squashed.xyz'
    string = smile(filename)
    vec = smiles_encoder(string)
    print(vec)


if __name__ == '__main__':
    test()
