# SELFIES

![versions](https://img.shields.io/pypi/pyversions/selfies.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)


SELFIES (SELF-referencIng Embedded Strings) is a general-purpose, 
sequence-based, robust representation of semantically constrained graphs. It
is based on a Chomsky type-2 grammar, augmented with two self-referencing 
functions. A main objective is to use SELFIES as direct input into machine 
learning models, in particular in generative models, for the generation of 
graphs with high semantical and syntactical validity (SELFIES has a validity 
of >99.99% even for random strings). The code presented here is a concrete 
application of SELFIES in chemistry, for the robust representation of molecule. 
See the paper at arXiv: https://arxiv.org/abs/1905.13741.


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

* More examples can be found in the ``examples/`` directory, including a 
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

By default, SELFIES is tested against a random 10,000 molecule subset of 
various datasets including:
 * [QM9](https://www.nature.com/articles/sdata201422), 
 * [ZINC](https://en.wikipedia.org/wiki/ZINC_database), 
 * [Non-fullerene acceptors for organic solar cells](https://www.sciencedirect.com/science/article/pii/S2542435117301307)
 * Tox21
 * PubChem MUV 


## License 

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
