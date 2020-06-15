## SELFIES

SELFIES (SELF-referencIng Embedded Strings) is a general-purpose, sequence-based,
robust representation of semantically constrained graphs. It is based on a Chomsky
type-2 grammar, augmented with two self-referencing functions. A main objective is
to use SELFIES as direct input into machine learning models, in particular in
generative models, for the generation of graphs with high semantical and syntactical
validity.

See the paper at arXiv: https://arxiv.org/abs/1905.13741

The code presented here is a concrete application of SELFIES in chemistry, for
the robust representation of molecule. 

SELFIES has a validity of >99.99% even for entire random strings. 

### Installation
You can install SELFIES via
```
pip install selfies
```

### Examples
Several examples can be seen in examples/selfies_example.py. Here is a simple encoding and decoding:

```python
from selfies import encoder, decoder, get_alphabet  
    
test_molecule1='CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1' # non-fullerene acceptors for organic solar cells
selfies1=encoder(test_molecule1)
smiles1=decoder(selfies1)

print('test_molecule1: '+test_molecule1+'\n')
print('selfies1: '+selfies1+'\n')
print('smiles1: '+smiles1+'\n')
print('equal: '+str(test_molecule1==smiles1)+'\n\n\n')

my_alphabet=get_alphabet() # contains all semantically valid SELFIES symbols.

```

- an example of SELFIES in a generative model can be seen in the directory 'VariationalAutoEncoder_with_SELFIES\'. There, SMILES datasets are automatically translated into SELFIES, and used for training of a variational autoencoder (VAE).

- One example used SELFIES in a genetic algorithm to achieve state-of-the-art performance for inverse design in this [ICLR2020 paper](https://arxiv.org/abs/1909.11655), with the [code here](https://github.com/aspuru-guzik-group/GA).

### Running Tests
SELFIES uses `pytest` as its testing framework. All tests can be found in the `tests/` directory.

You can run the test suite for SELFIES from your command line:

```bash
python setup.py test
```

These tests are necessary but not sufficient for the correctness of SELFIES. ToDo: More molecules should be tested, and the final comparison should be between the canonical input SMILES and the canonical output SMILES.

### Python version
fully tested with Python 3.7.1 on
- 134.000 molecules at QM9 database (https://www.nature.com/articles/sdata201422)
- 250.000 molecues from the ZINC database (https://en.wikipedia.org/wiki/ZINC_database)
- 72 million molecules from PubChem (https://pubchem.ncbi.nlm.nih.gov/)
- 50.000 molecules for organic solar cells (https://www.sciencedirect.com/science/article/pii/S2542435117301307)
- 1 million molecules from organic chemical reactions (https://pubs.rsc.org/en/content/articlehtml/2018/sc/c8sc02339e)

supported:
- Python 3.7.2, 3.7.1, 3.6.8, 3.6.7, 2.7.15



### Versions
#### 0.2.4 (01.10.2019):
       - added:
           -> functon get_alphabet() which returns a list of 29 selfies symbols whos arbitrary combination produce >99.99% valid molecules
       - bug fixes:
           -> fixed bug which happens when three rings start at one node, and two of them form a double ring
           -> enabled rings with sizes of up to 8000 SELFIES symbols
           -> bugfix for tiny ring to RDkit syntax conversion, spanning multiple branches
       - we thank Kevin Ryan (LeanAndMean@github), Theophile Gaudin and Andrew Brereton for suggestions and bug reports 

#### 0.2.2 (19.09.2019):
       - added:
           -> Enabled [C@],[C@H],[C@@],[C@@H],[H] to use in a semantic constrained way
       - we thank Andrew Brereton for suggestions and bug reports 


#### 0.2.1 (02.09.2019):
       - added:
           -> Decoder: added optional argument to restrict nitrogen to 3 bonds. decoder(...,N_restrict=False) to allow for more bonds;
                       standard: N_restrict=True
           -> Decoder: added optional argument make ring-function bi-local (i.e. confirms bond number at target).
                       decoder(...,bilocal_ring_function=False) to not allow bi-local ring function; standard:
                       bilocal_ring_function=True. The bi-local ring function will allow validity of >99.99% of random molecules
           -> Decoder: made double-bond ring RDKit syntax conform
           -> Decoder: added state X5 and X6 for having five and six bonds free
       - bug fixes:
            -> Decoder+Encoder: allowing for explicit brackets for organic atoms, for instance [I]
            -> Encoder: explicit single/double bond for non-canconical SMILES input issue fixed
            -> Decoder: bug fix for [Branch*] in state X1
       - we thank Benjamin Sanchez-Lengeling, Theophile Gaudin and Zhenpeng Yao for suggestions and bug reports 

#### 0.1.1 (04.06.2019): 
       - initial release 
