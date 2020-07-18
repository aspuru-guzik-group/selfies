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
application of SELFIES in chemistry, for the robust representation of
molecule. See the paper at arXiv (https://arxiv.org/abs/1905.13741)
by Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich, and
Alan Aspuru-Guzik.


## Installation
Use pip to install selfies.

```bash
pip install selfies
```

## Usage

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

## Documentation 

The documentation can be found on ReadTheDocs. Alternatively, it can be built 
from the ``docs/`` directory. 


## Tests
SELFIES uses `pytest` with `tox` as its testing framework.
All tests can be found in  the `tests/` directory. To run the test suite for 
SELFIES, install ``tox`` and run:  

```bash
tox
```

By default, SELFIES is tested against a random subset
(of size ``dataset_samples``) on various datasets:

 * 130K molecules from [QM9](https://www.nature.com/articles/sdata201422)
 * 250K molecules from [ZINC](https://en.wikipedia.org/wiki/ZINC_database), 
 * 50K molecules from [non-fullerene acceptors for organic solar cells](https://www.sciencedirect.com/science/article/pii/S2542435117301307)
 * 8K molecules from Tox21
 * 93K molecules from PubChem MUV
 
Other tests are random and repeated ``trials`` number of times. 
These can be specified as arguments 

```bash
tox -- --trials 100 --dataset_samples 100
```

where ``--trials=10000`` and ``--dataset_samples=10000`` by default. Note that 
if ``dataset_samples`` is negative or exceeds the length of the dataset, 
the whole dataset is used. 


## Credits

We thank Kevin Ryan (LeanAndMean@github), Theophile Gaudin, Andrew Brereton,
Benjamin Sanchez-Lengeling, and Zhenpeng Yao for their suggestions and 
bug reports. 

## License 

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
