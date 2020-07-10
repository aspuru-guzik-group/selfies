#
# Self-Referencing Embedded Strings (SELFIES):
#                  A 100% robust molecular string representation
#                  (https://arxiv.org/abs/1905.13741)
# by Mario Krenn, Florian Haese, AkshatKumar Nigam,
#      Pascal Friederich, Alán Aspuru-Guzik
#
# Demo of Rubustness of SMILES and SELFIES and DeepSMILES
# 
# Generates 1000 cases of 1, 2 or 3 mutations of small bio-molecule (MDMA).
# The alphabets are those that can translate the QM0 dataset (and have been used in all experiments in the paper)
#
#
# questions/remarks: mario.krenn@utoronto.ca or alan@aspuru.com
#
#                        11.03.2020
#
# Requirements: RDKit
#               selfies (pip install selfies)    
#               DeepSMILES (pip install --upgrade deepsmiles)            
#
#


from rdkit.Chem import MolFromSmiles
from rdkit import rdBase

from random import randint
from selfies import encoder, decoder  

import deepsmiles

rdBase.DisableLog('rdApp.error')


def IsCorrectSMILES(smiles):
    if len(smiles)==0:
        resMol=None
    else:
        try:
            resMol=MolFromSmiles(smiles, sanitize=True)
        except Exception:
            resMol=None

    if resMol==None:
        return 0
    else:
        return 1


def tokenize_selfies(selfies):
    location=selfies.find(']')
    all_tokens=[]
    while location>=0:
        all_tokens.append(selfies[0:location+1])
        selfies=selfies[location+1:]
        location=selfies.find(']')
        
    return all_tokens



def detokenize_selfies(selfies_list):
    selfies=''
    for ii in range(len(selfies_list)):
        selfies=selfies+selfies_list[ii]

    return selfies

mdma='CNC(C)CC1=CC=C2C(=C1)OCO2'
smiles_symbols='FONC()=#12345' # with this alphabet, the whole QM9 db can be translated (except of ions and stereochemistry)

print('\n\n\n')
print('SMILES: ',mdma,'\n')

num_repeat=1000
for c_num_of_mut in range(3):
    single_mut_err=0
    for c_muts in range(num_repeat):

        new_mdma=mdma
        for ii in range(c_num_of_mut+1):
            mol_idx=randint(0,len(new_mdma)-1)        
            symbol_idx=randint(0,len(smiles_symbols)-1)

            new_mdma=new_mdma[0:mol_idx]+smiles_symbols[symbol_idx]+new_mdma[mol_idx+1:]

        res_new=IsCorrectSMILES(new_mdma)
        if res_new==0:
            single_mut_err=single_mut_err+1

    print(c_num_of_mut+1, 'mutations with SMILES. Correct: ', num_repeat-single_mut_err, '/', num_repeat, '=', 1-single_mut_err/num_repeat)



# SELFIES code
mdma_selfies=encoder(mdma)
print('\n\n\n')
print('SELFIES: ',mdma_selfies,'\n')
mdma_selfies_tok=tokenize_selfies(mdma_selfies)
selfies_symbols=['[epsilon]','[Ring1]','[Ring2]','[Branch1_1]','[Branch1_2]','[Branch1_3]','[F]','[O]','[=O]','[N]','[=N]','[#N]','[C]','[=C]','[#C]']; #with this alphabet, the whole QM9 db can be translated (except of ions and stereochemistry)

num_repeat=1000
for c_num_of_mut in range(3):
    single_mut_err=0
    for c_muts in range(num_repeat):

        new_mdma=mdma_selfies_tok
        for ii in range(c_num_of_mut+1):
            mol_idx=randint(0,len(new_mdma)-1)        
            symbol_idx=randint(0,len(selfies_symbols)-1)

            new_mdma_str=detokenize_selfies(new_mdma[0:mol_idx])
            new_mdma_str=new_mdma_str+selfies_symbols[symbol_idx]
            new_mdma_str=new_mdma_str+detokenize_selfies(new_mdma[mol_idx+1:])
            new_mdma=tokenize_selfies(new_mdma_str)

        mutated_selfies=detokenize_selfies(new_mdma)
        mutated_smiles=decoder(mutated_selfies)
        res_new=IsCorrectSMILES(mutated_smiles)
        
        if res_new==0:
            single_mut_err=single_mut_err+1
            
        if c_muts>0 and c_muts%1000==0:
            print('Iteration: ', c_muts, '/', num_repeat)
        

    print(c_num_of_mut+1, 'mutations with SELFIES. Correct: ', num_repeat-single_mut_err, '/', num_repeat, '=', 1-single_mut_err/num_repeat)



# DeepSMILES code
deepsmiles_symbols='FONC)=#3456789' # with this alphabet, the whole QM9 db can be translated (except of ions and stereochemistry)
converter = deepsmiles.Converter(rings=True, branches=True)

mdma_deepsmiles=converter.encode(mdma)
print('\n\n\n')
print('DeepSMILES: ',mdma_deepsmiles,'\n')

num_repeat=1000
for c_num_of_mut in range(3):
    single_mut_err=0
    for c_muts in range(num_repeat):

        new_mdma=mdma_deepsmiles
        for ii in range(c_num_of_mut+1):
            mol_idx=randint(0,len(new_mdma)-1)        
            symbol_idx=randint(0,len(deepsmiles_symbols)-1)

            new_mdma=new_mdma[0:mol_idx]+deepsmiles_symbols[symbol_idx]+new_mdma[mol_idx+1:]

        try:
            mutated_smiles=converter.decode(new_mdma)
        except Exception:
            mutated_smiles='err'

        res_new=IsCorrectSMILES(mutated_smiles)
        if res_new==0:
            single_mut_err=single_mut_err+1

    print(c_num_of_mut+1, 'mutations with DeepSMILES. Correct: ', num_repeat-single_mut_err, '/', num_repeat, '=', 1-single_mut_err/num_repeat)
