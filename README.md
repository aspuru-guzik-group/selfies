## SELFIES

SELFIES (SELF-referencIng Embedded Strings) is a general-purpose, sequence-based,
robust representation of semantically constrained graphs. It is based on a Chomsky
type-2 grammar, augmented with two self-referencing functions. A main objective is
to use SELFIES as direct input into machine learning models, in particular
in generative models, for the generation of outputs with high validity.

The code presented here is a concrete application of SELIFES in chemistry, for
the robust representation of molecule. We show the encoding and decoding of three
molecules from various databases, and the generation of a new, random molecule
with high semantical and syntactical validity.

### Installation
You can install SELFIES via
```
pip install selfies
```

### Examples
Several examples can be seen in examples/selfies_example.py. Here is a simple encoding and decoding:

```python
from selfies import encoder, decoder
    
test_molecule1='CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1' # non-fullerene acceptors for organic solar cells
selfies1=encoder(test_molecule1)
smiles1=decoder(selfies1)

print('test_molecule1: '+test_molecule1+'\n')
print('selfies1: '+selfies1+'\n')
print('smiles1: '+smiles1+'\n')
print('equal: '+str(test_molecule1==smiles1)+'\n\n\n')
```

### Python version
fully tested with Python 3.7.1 on
- 134.000 molecules at QM9 database (https://www.nature.com/articles/sdata201422)
- 250.000 molecues from the ZINC database (https://en.wikipedia.org/wiki/ZINC_database)
- 72 million molecules from PubChem (https://pubchem.ncbi.nlm.nih.gov/)
- 50.000 molecules for organic solar cells (https://www.sciencedirect.com/science/article/pii/S2542435117301307)
- 1 million molecules from organic chemical reactions (https://pubs.rsc.org/en/content/articlehtml/2018/sc/c8sc02339e)

supported:
- Python 3.7.2, 3.7.1, 3.6.8, 3.6.7, 2.7.15
