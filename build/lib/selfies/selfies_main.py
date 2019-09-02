# =============================================================================
# 
# SELFIES: a robust representation of semantically constrained graphs with an example application in chemistry
#               v0.2.0, 02. September 2019
# by Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich, Alan Aspuru-Guzik
#
# SELFIES (SELF-referencIng Embedded Strings) is a general-purpose, sequence-based,
# robust representation of semantically constrained graphs. It is based on a Chomsky
# type-2 grammar, augmented with two self-referencing functions. A main objective is
# to use SELFIES as direct input into machine learning models, in particular
# in generative models, for the generation of outputs with high validity.
#
# The code presented here is a concrete application of SELIFES in chemistry, for
# the robust representation of molecule. We show the encoding and decoding of three
# molecules from various databases, and the generation of a new, random molecule
# with high semantical and syntactical validity.
#
# This file contains the encoder (SMILES -> SELFIES) and decoder (SELFIES -> SMILES),
# as well as an example for creating random SELFIES.
# 
#
# fully tested with Python 3.7.1 on
#     - 134.000 molecules at QM9 database (https://www.nature.com/articles/sdata201422)
#     - 250.000 molecues from the ZINC database (https://en.wikipedia.org/wiki/ZINC_database)
#     - 72 million molecules from PubChem (https://pubchem.ncbi.nlm.nih.gov/)
#     - 50.000 molecules for organic solar cells (https://www.sciencedirect.com/science/article/pii/S2542435117301307)
#     - 1 million molecules from organic chemical reactions (https://pubs.rsc.org/en/content/articlehtml/2018/sc/c8sc02339e)
#
# supported:
# - Python 3.7.2
# - Python 3.7.1
# - Python 3.6.8
# - Python 3.6.7
# - Python 2.7.15
#
#
# versions:
# 0.2.0 (02.09.2019):
#       - added:
#           -> Decoder: added optional argument to restrict nitrogen to 3 bonds. decoder(...,N_restrict=False) to allow for more bonds;
#                       standard: N_restrict=True
#           -> Decoder: added optional argument make ring-function bi-local (i.e. confirms bond number at target).
#                       decoder(...,bilocal_ring_function=False) to not allow bi-local ring function; standard:
#                       bilocal_ring_function=True. The bi-local ring function will allow validity of >99.99% of random molecules
#           -> Decoder: made double-bond ring RDKit syntax conform
#           -> Decoder: added state X5 and X6 for having five and six bonds free
#       - bug fixes:
#            -> Decoder+Encoder: allowing for explicit brackets for organic atoms, for instance [I]
#            -> Encoder: explicit single/double bond for non-canconical SMILES input issue fixed
#            -> Decoder: bug fix for [Branch*] in state X1
#       - we thank Benjamin Sanchez-Lengeling, Theophile Gaudin and Zhenpeng Yao for suggestions and bug reports 
#
# 0.1.1 (04.06.2019): 
#       - initial release    
#
#
# For comments, bug reports or feature ideas, please send an email to
# mario.krenn@utoronto.ca and alan@aspuru.com
# =============================================================================

from random import randint
from selfies_fcts import encoder, decoder

# Now we encode three molecules from SMILES -> SELFIES, and decode them from SELFIES -> SMILES
test_molecule1='CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1' # non-fullerene acceptors for organic solar cells
selfies1=encoder(test_molecule1)
smiles1=decoder(selfies1)
print('test_molecule1: '+test_molecule1+'\n')
print('selfies1: '+selfies1+'\n')
print('smiles1: '+smiles1+'\n')
print('equal: '+str(test_molecule1==smiles1)+'\n\n\n')

test_molecule2='CC(C)c1noc(-c2cc[nH+]c(N3CCN(C(=O)[C@H]4C[C@H]4C)CC3)c2)n1' # from ZINC database
selfies2=encoder(test_molecule2)
smiles2=decoder(selfies2)
print('test_molecule2: '+test_molecule2+'\n')
print('selfies2: '+selfies2+'\n')
print('smiles2: '+smiles2+'\n')
print('equal: '+str(test_molecule2==smiles2)+'\n\n\n')


test_molecule3='CCOC(=O)C1(C(=O)OCC)C23c4c5c6c7c8c4-c4c2c2c9c%10c4C4%11c%12c-%10c%10c%13c%14c%15c%16c%17c%18c%19c%20c%21c%22c%23c%24c(c-7c(c7c%12c%13c(c7%24)c(c%19%23)c%18%14)C84C%11(C(=O)OCC)C(=O)OCC)C%224C(C(=O)OCC)(C(=O)OCC)C64c4c-5c5c6c(c4-%21)C%204C(C(=O)OCC)(C(=O)OCC)C%174c4c-6c(c-2c(c4-%16)C92C(C(=O)OCC)(C(=O)OCC)C%10%152)C513' # from PubChem
selfies3=encoder(test_molecule3)
smiles3=decoder(selfies3)
print('test_molecule1: '+test_molecule3+'\n')
print('selfies1: '+selfies3+'\n')
print('smiles1: '+smiles3+'\n')
print('equal: '+str(test_molecule3==smiles3)+'\n\n\n')


#Create a random Molecule
my_alphabet=['[Branch1_1]','[Branch1_2]','[Branch1_3]','[Ring1]','[O]','[=O]','[N]','[=N]','[C]','[=C]','[#C]','[S]','[=S]','[P]','[F]']; # this is a very small alphabet from which the random selfies are generated
                                                                                                                                          # This alphabet can be extended with additional elements. For example, see the list start_alphabet in the function smiles_to_selfies.
                                                                                                                                          # Also when you run the three test-molecules above, you see the brackets that are used, and can use some of them.
                                                                                                                                          
len_of_molecule=30 # Number of selfies symbols of the random string. The final SMILES string will not necessarily be of the same size, because some elements of this alphabet stop the derivation (such as Flour, as it can form only a single bond)
                  
rnd_selfies=''
for ii in range(len_of_molecule):
    rnd_selfies+=my_alphabet[randint(0,len(my_alphabet)-1)]

smiles4=decoder(rnd_selfies)
print('From random selfies: '+smiles4+'\n')

