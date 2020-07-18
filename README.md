# SELFIES

![versions](https://img.shields.io/pypi/pyversions/selfies.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)


SELFIES (SELF-referencIng Embedded Strings) is a general-purpose, 
sequence-based, robust representation of semantically constrained graphs. It
is based on a Chomsky type-2 grammar, augmented with two self-referencing 
functions. A main objective is to use SELFIES as direct input into machine 
learning models, in particular in generative models, for the generation of 
graphs with high semantical and syntactical validity (SELFIES has a validity 
of >99.99% even for entire random strings). 

The code presented here is a concrete application of SELFIES in chemistry, for
the robust representation of molecule. See the paper at
[arXiv](https://arxiv.org/abs/1905.13741).


## Installation
Use pip to install selfies.

```bash
pip install selfies
```

## Code Examples

```python
import selfies as sf
    
benzene = "C1=CC=CC=C1"
encoded_selfies = sf.encoder(benzene)  # SMILES --> SEFLIES
decoded_smiles = sf.decoder(encoded_selfies)  # SELFIES --> SMILES
```

* More examples can be seen in the ``examples/`` directory, including a 
variational autoencoder that runs on the SELFIES language.
* This [ICLR2020 paper](https://arxiv.org/abs/1909.11655) used SELFIES in a
genetic algorithm to achieve state-of-the-art performance for inverse design, 
with the [code here](https://github.com/aspuru-guzik-group/GA).

## Tests
SELFIES uses `pytest` as its testing framework. All tests can be found in 
the `tests/` directory.

Many of the SELFIES tests use [RDKit](https://www.rdkit.org/), which can 
be installed using Conda. To run the test suite for SELFIES, create a Conda
environment with RDKit installed, and run from your command line:  

```bash
python setup.py test
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


## License 

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
