# =============================================================================
# 
# SELFIES: a robust representation of semantically constrained graphs with an example application in chemistry
# Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich, Alan Aspuru-Guzik
# 
# This file contains the encoder (SMILES -> SELFIES) and decoder (SELFIES -> SMILES),
# as well as an example for creating random SELFIES.
# 
# For comments, bug reports or feature ideas, please send an email to
# mario.krenn@utoronto.ca and alan@aspuru.com
#
# v1.0, 01.06.2019
#
# 
# =============================================================================

from random import randint
from selfies_fcts import make_brackets_around_atoms, reconfigure_smiles_numbers1, reconfigure_smiles_numbers2, smiles_to_selfies, selfies_to_smiles, insert_rings_to_smiles

def selfies_encoder(smiles,PrintErrorMessage=True):
    try:
        preselfies1=make_brackets_around_atoms(smiles)         # Itemize Atoms
        preselfies2=reconfigure_smiles_numbers1(preselfies1)  # Unique Ringsymbols
        preselfies3=reconfigure_smiles_numbers2(preselfies2)  # Rings as relative distance
        selfies=smiles_to_selfies(preselfies3)                # Create selfies
    except ValueError as err:
        if PrintErrorMessage:
            print(err)
            print('Could not encode smiles string. Please contact authors.')
        return -1
        
    return selfies
    
def selfies_decoder(selfies,PrintErrorMessage=True):
    smiles=-1
    if selfies!=-1:
        try:
            presmiles1=selfies_to_smiles(selfies)    # Runs Grammar Rules
            smiles=insert_rings_to_smiles(presmiles1) # Inserts Rings
        except  ValueError as err:
            if PrintErrorMessage:
                print(err)            
                print('Could not decode selfies string. Please contact authors.')
            return -1
    
    return smiles
    

test_molecule1='CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1' # non-fullerene acceptors for organic solar cells
selfies1=selfies_encoder(test_molecule1)
smiles1=selfies_decoder(selfies1)
print('test_molecule1: '+test_molecule1+'\n')
print('selfies1: '+selfies1+'\n')
print('smiles1: '+smiles1+'\n')
print('equal: '+str(test_molecule1==smiles1)+'\n\n\n')

test_molecule2='CC(C)c1noc(-c2cc[nH+]c(N3CCN(C(=O)[C@H]4C[C@H]4C)CC3)c2)n1' # from ZINC database
selfies2=selfies_encoder(test_molecule2)
smiles2=selfies_decoder(selfies2)
print('test_molecule2: '+test_molecule2+'\n')
print('selfies2: '+selfies2+'\n')
print('smiles2: '+smiles2+'\n')
print('equal: '+str(test_molecule2==smiles2)+'\n\n\n')


test_molecule3='CCOC(=O)C1(C(=O)OCC)C23c4c5c6c7c8c4-c4c2c2c9c%10c4C4%11c%12c-%10c%10c%13c%14c%15c%16c%17c%18c%19c%20c%21c%22c%23c%24c(c-7c(c7c%12c%13c(c7%24)c(c%19%23)c%18%14)C84C%11(C(=O)OCC)C(=O)OCC)C%224C(C(=O)OCC)(C(=O)OCC)C64c4c-5c5c6c(c4-%21)C%204C(C(=O)OCC)(C(=O)OCC)C%174c4c-6c(c-2c(c4-%16)C92C(C(=O)OCC)(C(=O)OCC)C%10%152)C513' # from PubChem
selfies3=selfies_encoder(test_molecule3)
smiles3=selfies_decoder(selfies3)
print('test_molecule1: '+test_molecule3+'\n')
print('selfies1: '+selfies3+'\n')
print('smiles1: '+smiles3+'\n')
print('equal: '+str(test_molecule3==smiles3)+'\n\n\n')


#Create a random Molecule
my_alphabet=['[Branch1_1]','[Branch1_2]','[Branch1_3]','[Ring1]','[O]','[=O]','[N]','[=N]','[C]','[=C]','[#C]','[S]','[=S]','[P]','[F]']; # this is a very small alphabet from which the random selfies are generated
                                                                                                                                          # This alphabet can be extended with additional elements. For example, see the list start_alphabet in the function smiles_to_selfies.
                                                                                                                                          # Also when you run the three test-molecules above, you see the brackets that are used, and can use some of them.
                                                                                                                                          
len_of_molecule=15 # Number of selfies symbols of the random string. The final SMILES string will not necessarily be of the same size, because some elements of this alphabet stop the derivation (such as Flour, as it can form only a single bond)
                  
rnd_selfies=''
for ii in range(len_of_molecule):
    rnd_selfies+=my_alphabet[randint(0,len(my_alphabet)-1)]

smiles4=selfies_decoder(rnd_selfies)
print('From random selfies: '+smiles4+'\n')

