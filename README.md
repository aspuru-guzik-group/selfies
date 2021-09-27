# SELFIES

[![GitHub release](https://img.shields.io/github/release/aspuru-guzik-group/selfies.svg)](https://GitHub.com/aspuru-guzik-group/selfies/releases/)
![versions](https://img.shields.io/pypi/pyversions/selfies.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/aspuru-guzik-group/selfies/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/aspuru-guzik-group/selfies.svg)](https://GitHub.com/aspuru-guzik-group/selfies/issues/)
[![Documentation Status](https://readthedocs.org/projects/selfies/badge/?version=latest)](http://selfies.readthedocs.io/?badge=latest)
[![GitHub contributors](https://img.shields.io/github/contributors/aspuru-guzik-group/selfies.svg)](https://GitHub.com/aspuru-guzik-group/selfies/graphs/contributors/)


**Self-Referencing Embedded Strings (SELFIES): A 100% robust molecular string representation**\
_Mario Krenn, Florian Haese, AkshatKumar Nigam, Pascal Friederich, Alan Aspuru-Guzik_\
[*Machine Learning: Science and Technology* **1**, 045024 (2020)](https://iopscience.iop.org/article/10.1088/2632-2153/aba947), [extensive blog post January 2021](https://aspuru.substack.com/p/molecular-graph-representations-and).\
[Talk on youtube about SELFIES](https://www.youtube.com/watch?v=CaIyUmfGXDk).\
[Blog explaining SELFIES in Japanese language](https://blacktanktop.hatenablog.com/entry/2021/08/12/115613)\
Major contributors since v1.0.0: _[Alston Lo](https://github.com/alstonlo) and [Seyone Chithrananda](https://github.com/seyonechithrananda)_\
Chemistry Advisor: [Robert Pollice](https://scholar.google.at/citations?user=JR2N3JIAAAAJ)

---

A main objective is to use SELFIES as direct input into machine learning models,
in particular in generative models, for the generation of molecular graphs
which are syntactically and semantically valid.

<p align="center">
   <img src="https://github.com/aspuru-guzik-group/selfies/blob/master/examples/VAE_LS_Validity.png" alt="SELFIES validity in a VAE latent space" width="666px">
</p>

## Installation
Use pip to install ``selfies``.

```bash
pip install selfies
```

To check if the correct version of ``selfies`` is installed, use
the following pip command. 

```bash
pip show selfies
```

To upgrade to the latest release of ``selfies`` if you are using an 
older version, use the following pip command. Please see the 
[CHANGELOG](https://github.com/aspuru-guzik-group/selfies/blob/master/CHANGELOG.md) 
to review the changes between versions of `selfies`, before upgrading: 

```bash
pip install selfies --upgrade 
```

## Documentation

The documentation can be found on
[ReadTheDocs](https://selfies.readthedocs.io/en/latest/).
Alternatively, it can be built from the ``docs/`` directory.

## Usage

### Standard Functions

The ``selfies`` library has eight standard functions:

| Function | Description |
| -------- | ----------- |
| ``selfies.encoder`` | Translates a SMILES into an equivalent SELFIES. |
| ``selfies.decoder`` | Translates a SELFIES into an equivalent SMILES. |
| ``selfies.len_selfies`` | Returns the (symbol) length of a SELFIES.  |
| ``selfies.split_selfies`` | Splits a SELFIES into its symbols. |
| ``selfies.get_alphabet_from_selfies`` | Builds an alphabet of SELFIES symbols from an iterable of SELFIES. |
| ``selfies.get_semantic_robust_alphabet`` | Returns a subset of all SELFIES symbols that are semantically constrained. |
| ``selfies.selfies_to_encoding`` | Converts a SELFIES into a label and/or one-hot encoding. |
| ``selfies.encoding_to_selfies`` | Converts a label or one-hot encoding into a SELFIES. |

Please read the documentation for more detailed descriptions of these
functions, and to view the advanced functions, which allow users to
customize the SELFIES language.

### Examples

#### Translation between SELFIES and SMILES representations:

```python
import selfies as sf

benzene = "c1ccccc1"

# SMILES -> SELFIES -> SMILES translation
try:
    benzene_sf = sf.encoder(benzene)  # [C][=C][C][=C][C][=C][Ring1][=Branch1]
    benzene_smi = sf.decoder(benzene_sf)  # C1=CC=CC=C1
except sf.EncoderError:
    pass  # sf.encoder error!
except sf.DecoderError:
    pass  # sf.decoder error!

len_benzene = sf.len_selfies(benzene_sf)  # 8

symbols_benzene = list(sf.split_selfies(benzene_sf))
# ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[=Branch1]']
```

#### Customizing SELFIES


#### Integer and one-hot encoding SELFIES:
In this example we first build an alphabet
from a dataset of SELFIES, and then convert a SELFIES into a
padded, label-encoded representation. Note that we use the
``[nop]`` ([no operation](https://en.wikipedia.org/wiki/NOP_(code) ))
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

# SELFIES to label encode
dimethyl_ether = dataset[0]  # '[C][O][C]'

# [1, 3, 1, 4, 4]
print(sf.selfies_to_encoding(dimethyl_ether,
                             vocab_stoi=symbol_to_idx,
                             pad_to_len=pad_to_len,
                             enc_type='label'))
                             
# [[0, 1, 0, 0, 0], [0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 1]]
print(sf.selfies_to_encoding(dimethyl_ether,
                             vocab_stoi=symbol_to_idx,
                             pad_to_len=pad_to_len,
                             enc_type='one_hot'))
```

### More Examples

* More examples can be found in the ``examples/`` directory, including a
[variational autoencoder that runs on the SELFIES](https://github.com/aspuru-guzik-group/selfies/tree/master/examples/vae_example) language.
* This [ICLR2020 paper](https://arxiv.org/abs/1909.11655) used SELFIES in a
genetic algorithm to achieve state-of-the-art performance for inverse design,
with the [code here](https://github.com/aspuru-guzik-group/GA).
* SELFIES allows for [highly efficient exploration and interpolation of the chemical space](https://chemrxiv.org/articles/preprint/Beyond_Generative_Models_Superfast_Traversal_Optimization_Novelty_Exploration_and_Discovery_STONED_Algorithm_for_Molecules_using_SELFIES/13383266), with a [deterministic algorithms, see code](https://github.com/aspuru-guzik-group/stoned-selfies).
* We use SELFIES for [Deep Molecular dreaming](https://arxiv.org/abs/2012.09712), a new generative model inspired by interpretable neural networks in computational vision. See the [code of PASITHEA here](https://github.com/aspuru-guzik-group/Pasithea).
* Kohulan Rajan, Achim Zielesny, Christoph Steinbeck show in two papers that SELFIES outperforms other representations in [img2string](https://link.springer.com/article/10.1186/s13321-020-00469-w) and [string2string](https://chemrxiv.org/articles/preprint/STOUT_SMILES_to_IUPAC_Names_Using_Neural_Machine_Translation/13469202/1) translation tasks, see the codes of [DECIMER](https://github.com/Kohulan/DECIMER-Image-to-SMILES) and [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator). 
* An improvement to the old genetic algorithm, the authors have also released [JANUS](https://arxiv.org/abs/2106.04011), which allows for more efficient optimization in the chemical space. JANUS makes use of [STONED-SELFIES](https://pubs.rsc.org/en/content/articlepdf/2021/sc/d1sc00231g) and a neural network for efficient sampling. 

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
 * 250K molecules from [ZINC](https://en.wikipedia.org/wiki/ZINC_database)
 * 50K molecules from [non-fullerene acceptors for organic solar cells](https://www.sciencedirect.com/science/article/pii/S2542435117301307)
 * 160K+ molecules from various [MoleculeNet](http://moleculenet.ai/datasets-1) datasets
 * 36M+ molecules from the [eMolecules Database](https://www.emolecules.com/info/products-data-downloads.html).
   Due to its large size, this dataset is not included on the repository. To run tests 
   on it, please download the dataset in the ``tests/test_sets`` directory 
   and run the ``tests/run_on_large_dataset.py`` script. 

Other tests are random and repeated ``trials`` number of times.
These can be specified as arguments

```bash
tox -- --trials 100 --dataset_samples 100
```

where ``--trials=100000`` and ``--dataset_samples=100000`` by default. Note that
if ``dataset_samples`` is negative or exceeds the length of the dataset,
the whole dataset is used.

## Version History
See [CHANGELOG](https://github.com/aspuru-guzik-group/selfies/blob/master/CHANGELOG.md).

## Credits

We thank Jacques Boitreaud, Andrew Brereton, Matthew Carbone (x94carbone), Nathan Frey (ncfrey), Theophile Gaudin,
HelloJocelynLu, Hyunmin Kim (hmkim), Minjie Li, Vincent Mallet, Alexander Minidis (DocMinus), Kohulan Rajan (Kohulan),
Kevin Ryan (LeanAndMean), Benjamin Sanchez-Lengeling, and Zhenpeng Yao for their suggestions and bug reports,
and Robert Pollice for chemistry advices.

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
