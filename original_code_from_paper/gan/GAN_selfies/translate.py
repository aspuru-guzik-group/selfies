#alphabets = ['I', 'N', 'D', 'K', 'L', 'B', 'C', 'J', 'G', 'A', 'H', 'E', 'F', ' ']
#smile = IncludeRingsForSMILES(GrammarPlusToSMILES(smile,'X0'))

import time
from adjusted_selfies_fcts import encoder, decoder
import pandas as pd
import numpy as np

from rdkit.Chem import MolFromSmiles, MolToSmiles


def is_correct_smiles(smiles):    
    """
    Using RDKit to calculate whether molecule is syntactically and semantically valid.
    """
    try:
        res_molecule=MolFromSmiles(smiles, sanitize=True)
    except Exception:
        res_molecule=None
                
    if res_molecule==None:
        return 0
    else:
        return 1

def old_vs_new_decoder(old_selfies):
    translate_alphabet={'A': '[epsilon]', 'B': '[F]', 'C': '[=O]', 'D': '[#N]', 'E': '[O]', 'F': '[N]', 'G': '[=N]', 'H': '[C]', 'I': '[=C]', 'J': '[#C]', 'K': '[Branch1_1]', 'L': '[Branch1_3]', 'M': '[Branch1_2]', 'N': '[Ring1]'}
    new_selfie=''
    for jj in range(len(old_selfies)):
        new_selfie+=translate_alphabet[old_selfies[jj]]
    translated_mol=decoder(new_selfie)
    return translated_mol



df = pd.read_csv('0SelectedSMILES_QM9.txt')
smiles_list = np.asanyarray(df.smiles)
    
df = pd.read_csv('2RGSMILES_QM9.txt')
old_selfies_list = np.asanyarray(df.smiles)    

for ii in range(len(smiles_list)):
    new_smiles=old_vs_new_decoder(old_selfies_list[ii])
    isSame=(new_smiles==smiles_list[ii])
    if not isSame:
        print(str(ii)+': '+new_smiles+', '+smiles_list[ii])
        time.sleep(0.5)
        
    if ii%1000==0 and ii>0:
        print(ii,': ',new_smiles)