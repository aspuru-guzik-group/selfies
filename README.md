# SELFIES

![versions](https://img.shields.io/pypi/pyversions/selfies.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)


SELFIES (SELF-referencIng Embedded Strings) is a 100% robust molecular
string representation.

A main objective is to use SELFIES as direct input into machine learning
models, in particular in generative models, for the generation of molecular
graphs which are syntactically and semantically valid.

See the paper by Mario Krenn, Florian Haese, AkshatKumar Nigam,
Pascal Friederich, and Alan Aspuru-Guzik at
arXiv (https://arxiv.org/abs/1905.13741).


## Installation
Use pip to install ``selfies``.

```bash
pip install selfies
```

## Usage

### Standard Functions

The ``selfies`` library has six standard functions:

| Function | Description |
| -------- | ----------- |
| ``selfies.encoder`` | Translates a SMILES into an equivalent SELFIES. |
| ``selfies.decoder`` | Translates a SELFIES into an equivalent SMILES. |
| ``selfies.len_selfies`` | Returns the (symbol) length of a SELFIES.  |
| ``selfies.split_selfies`` | Splits a SELFIES into its symbols. |
| ``selfies.get_alphabet_from_selfies`` | Builds an alphabet of SELFIES symbols from an iterable of SELFIES. |
| ``selfies.get_semantic_robust_alphabet`` | Returns a subset of all SELFIES symbols that are semantically constrained. |

Please read the documentation for more detailed descriptions of these
functions, and to view the advanced functions, which allow users to
customize the SELFIES language.

### Examples

#### Translation between SELFIES and SMILES representations:

```python
import selfies as sf

benzene = "c1ccccc1"

# SMILES --> SELFIES translation
encoded_selfies = sf.encoder(benzene)  # '[C][=C][C][=C][C][=C][Ring1][Branch1_2]'

# SELFIES --> SMILES translation
decoded_smiles = sf.decoder(encoded_selfies)  # 'C1=CC=CC=C1'

len_benzene = sf.len_selfies(encoded_selfies)  # 8

symbols_benzene = list(sf.split_selfies(encoded_selfies))
# ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[Branch1_2]']
```

#### Integer encoding SELFIES:
In this example we first build an alphabet
from a dataset of SELFIES, and then convert a SELFIES into a
padded, integer-encoded representation. Note that we use the
``'[nop]'`` ([no operation](https://en.wikipedia.org/wiki/NOP_(code) ))
symbol to pad our SELFIES, which is a special SELFIES symbol that is always
ignored and skipped over by ``selfies.decoder``, making it a useful
padding character.

```python
import selfies as sf

dataset = ['[C][O][C]', '[F][C][F]', '[O][=O]', '[C][C][O][C][C]']
alphabet = sf.get_alphabet_from_selfies(dataset)
alphabet.add('[nop]')  # '[nop]' is a special padding symbol
alphabet = list(sorted(alphabet))
print(alphabet)  # ['[=O]', '[C]', '[F]', '[O]', '[nop]']

pad_to_len = max(sf.len_selfies(s) for s in dataset)  # 5
symbol_to_idx = {s: i for i, s in enumerate(alphabet)}

# SELFIES to integer encode
dimethyl_ether = dataset[0]  # '[C][O][C]'

# pad the SELFIES
dimethyl_ether += '[nop]' * (pad_to_len - sf.len_selfies(dimethyl_ether))

# integer encode the SELFIES
int_encoded = []
for symbol in sf.split_selfies(dimethyl_ether):
    int_encoded.append(symbol_to_idx[symbol])

print(int_encoded)  # [1, 3, 1, 4, 4]
```

### More Examples

* More examples can be found in the ``examples/`` directory, including a
variational autoencoder that runs on the SELFIES language.
* This [ICLR2020 paper](https://arxiv.org/abs/1909.11655) used SELFIES in a
genetic algorithm to achieve state-of-the-art performance for inverse design,
with the [code here](https://github.com/aspuru-guzik-group/GA).

## Documentation

The documentation can be found on
[ReadTheDocs](https://selfies-mirror.readthedocs.io/en/latest/?).
Alternatively, it can be built from the ``docs/`` directory.


## Tests
SELFIES uses `pytest` with `tox` as its testing framework.
All tests can be found in  the `tests/` directory. To run the test suite for
SELFIES, install ``tox`` and run:  

```bash
tox
```

By default, SELFIES is tested against a random subset
(of size ``dataset_samples=100000``) on various datasets:

 * 130K molecules from [QM9](https://www.nature.com/articles/sdata201422)
 * 250K molecules from [ZINC](https://en.wikipedia.org/wiki/ZINC_database),
 * 50K molecules from [non-fullerene acceptors for organic solar cells](https://www.sciencedirect.com/science/article/pii/S2542435117301307)
 * 8K molecules from [Tox21](http://moleculenet.ai/datasets-1) in MoleculeNet
 * 93K molecules from PubChem [MUV](http://moleculenet.ai/datasets-1) in MoleculeNet
 * 6M molecules from the [eMolecules Plus Database](https://www.emolecules.com/info/plus/download-database) (∗) 

Other tests are random and repeated ``trials`` number of times.
These can be specified as arguments

```bash
tox -- --trials 100 --dataset_samples 100
```

where ``--trials=100000`` and ``--dataset_samples=100000`` by default. Note that
if ``dataset_samples`` is negative or exceeds the length of the dataset,
the whole dataset is used.

---

(∗) The eMolecules dataset is not included in the GitHub repository, due to its size.
For developers interested in testing on this dataset, one can download it
and then split the file into chunks of 1 million SMILES
using the following command in Windows:

```
split -l 1000000 -d --additional-suffix=.txt version.txt split
```

On Mac, the following commands can be used:
```
brew install coreutils
gsplit -l 1000000 -d --additional-suffix=.txt version.txt split
```

Once the dataset is downloaded and split in the `tests/test_sets/` directory,
`tests/test_on_emolecule.py` can be run to test the eMolecules dataset. 


## Credits

We thank Kevin Ryan (LeanAndMean@github), Theophile Gaudin, Andrew Brereton,
Benjamin Sanchez-Lengeling, and Zhenpeng Yao for their suggestions and
bug reports.

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
