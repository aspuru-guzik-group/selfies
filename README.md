# SELFIES

[![GitHub release](https://img.shields.io/github/release/aspuru-guzik-group/selfies.svg)](https://GitHub.com/aspuru-guzik-group/selfies/releases/)
![versions](https://img.shields.io/pypi/pyversions/selfies.svg)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/aspuru-guzik-group/selfies/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/aspuru-guzik-group/selfies.svg)](https://GitHub.com/aspuru-guzik-group/selfies/issues/)
[![Documentation Status](https://readthedocs.org/projects/selfiesv2/badge/?version=latest)](http://selfiesv2.readthedocs.io/?badge=latest)
[![GitHub contributors](https://img.shields.io/github/contributors/aspuru-guzik-group/selfies.svg)](https://GitHub.com/aspuru-guzik-group/selfies/graphs/contributors/)


**Self-Referencing Embedded Strings (SELFIES): A 100% robust molecular string representation**\
_Mario Krenn, Florian Haese, AkshatKumar Nigam, Pascal Friederich, Alan Aspuru-Guzik_\
[*Machine Learning: Science and Technology* **1**, 045024 (2020)](https://iopscience.iop.org/article/10.1088/2632-2153/aba947), [extensive blog post January 2021](https://aspuru.substack.com/p/molecular-graph-representations-and).\
[Talk on youtube about SELFIES](https://www.youtube.com/watch?v=CaIyUmfGXDk).\
[A community paper with 31 authors on SELFIES and the future of molecular string representations](https://arxiv.org/abs/2204.00056).\
[Blog explaining SELFIES in Japanese language](https://blacktanktop.hatenablog.com/entry/2021/08/12/115613)\
Major contributors of v1.0.n: _[Alston Lo](https://github.com/alstonlo) and [Seyone Chithrananda](https://github.com/seyonechithrananda)_\
Main developer of v2.0.0: _[Alston Lo](https://github.com/alstonlo)_\
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


## Usage

### Overview

Please refer to the [documentation](https://selfiesv2.readthedocs.io/en/latest/),
which contains a thorough tutorial  for getting started with ``selfies``
and detailed descriptions of the functions
that ``selfies`` provides. We summarize some key functions below.

| Function                              | Description                                                       |
| ------------------------------------- | ----------------------------------------------------------------- |
| ``selfies.encoder``                   | Translates a SMILES string into its corresponding SELFIES string. |
| ``selfies.decoder``                   | Translates a SELFIES string into its corresponding SMILES string. |
| ``selfies.set_semantic_constraints``  | Configures the semantic constraints that ``selfies`` operates on. |
| ``selfies.len_selfies``               | Returns the number of symbols in a SELFIES string.                |
| ``selfies.split_selfies``             | Tokenizes a SELFIES string into its individual symbols.           |
| ``selfies.get_alphabet_from_selfies`` | Constructs an alphabet from an iterable of SELFIES strings.       |
| ``selfies.selfies_to_encoding``       | Converts a SELFIES string into its label and/or one-hot encoding. |
| ``selfies.encoding_to_selfies``       | Converts a label or one-hot encoding into a SELFIES string.       |


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

#### Very simple creation of random valid molecules:
A key property of SELFIES is the possibility to create valid random molecules in a very simple way -- inspired by a tweet by [Rajarshi Guha](https://twitter.com/rguha/status/1543601839983284224):

```python
import selfies as sf
import random

alphabet=sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols
rnd_selfies=''.join(random.sample(list(alphabet), 9))
rnd_smiles=sf.decoder(rnd_selfies)
print(rnd_smiles)
```
These simple lines gives crazy molecules, but all are valid. Can be used as a start for more advanced filtering techniques or for machine learning models.

#### Integer and one-hot encoding SELFIES:

In this example, we first build an alphabet from a dataset of SELFIES strings,
and then convert a SELFIES string into its padded encoding. Note that we use the
``[nop]`` ([no operation](https://en.wikipedia.org/wiki/NOP_(code) ))
symbol to pad our SELFIES, which is a special SELFIES symbol that is always
ignored and skipped over by ``selfies.decoder``, making it a useful
padding character.

```python
import selfies as sf

dataset = ["[C][O][C]", "[F][C][F]", "[O][=O]", "[C][C][O][C][C]"]
alphabet = sf.get_alphabet_from_selfies(dataset)
alphabet.add("[nop]")  # [nop] is a special padding symbol
alphabet = list(sorted(alphabet))  # ['[=O]', '[C]', '[F]', '[O]', '[nop]']

pad_to_len = max(sf.len_selfies(s) for s in dataset)  # 5
symbol_to_idx = {s: i for i, s in enumerate(alphabet)}

dimethyl_ether = dataset[0]  # [C][O][C]

label, one_hot = sf.selfies_to_encoding(
   selfies=dimethyl_ether,
   vocab_stoi=symbol_to_idx,
   pad_to_len=pad_to_len,
   enc_type="both"
)
# label = [1, 3, 1, 4, 4]
# one_hot = [[0, 1, 0, 0, 0], [0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 1]]
```

#### Customizing SELFIES:

In this example, we relax the semantic constraints of ``selfies`` to allow
for hypervalences (caution: hypervalence rules are much less understood
than octet rules. Some molecules containing hypervalences are important,
but generally, it is not known which molecules are stable and reasonable).

```python
import selfies as sf

hypervalent_sf = sf.encoder('O=I(O)(O)(O)(O)O', strict=False)  # orthoperiodic acid
standard_derived_smi = sf.decoder(hypervalent_sf)
# OI (the default constraints for I allows for only 1 bond)

sf.set_semantic_constraints("hypervalent")
relaxed_derived_smi = sf.decoder(hypervalent_sf)
# O=I(O)(O)(O)(O)O (the hypervalent constraints for I allows for 7 bonds)
```

#### Explaining Translation:

You can get an "attribution" list that traces the connection between input and output tokens. For example let's see which tokens in the SELFIES string ``[C][N][C][Branch1][C][P][C][C][Ring1][=Branch1]`` are responsible for the output SMILES tokens.

```python
selfies = "[C][N][C][Branch1][C][P][C][C][Ring1][=Branch1]"
smiles, attr = sf.decoder(
    selfies, attribute=True)
print('SELFIES', selfies)
print('SMILES', smiles)
print('Attribution:')
for smiles_token in attr:
    print(smiles_token)

# output
SELFIES [C][N][C][Branch1][C][P][C][C][Ring1][=Branch1]
SMILES C1NC(P)CC1
Attribution:
AttributionMap(index=0, token='C', attribution=[Attribution(index=0, token='[C]')])
AttributionMap(index=2, token='N', attribution=[Attribution(index=1, token='[N]')])
AttributionMap(index=3, token='C', attribution=[Attribution(index=2, token='[C]')])
AttributionMap(index=5, token='P', attribution=[Attribution(index=3, token='[Branch1]'), Attribution(index=5, token='[P]')])
AttributionMap(index=7, token='C', attribution=[Attribution(index=6, token='[C]')])
AttributionMap(index=8, token='C', attribution=[Attribution(index=7, token='[C]')])
```

``attr`` is a list of `AttributionMap`s containing the output token, its index, and input tokens that led to it. For example, the ``P`` appearing in the output SMILES at that location is a result of both the ``[Branch1]`` token at position 3 and the ``[P]`` token at index 5. This works for both encoding and decoding. For finer control of tracking the translation (like tracking rings), you can access attributions in the underlying molecular graph with ``get_attribution``.

### More Usages and Examples

* More examples can be found in the ``examples/`` directory, including a
[variational autoencoder that runs on the SELFIES](https://github.com/aspuru-guzik-group/selfies/tree/master/examples/vae_example) language.
* This [ICLR2020 paper](https://arxiv.org/abs/1909.11655) used SELFIES in a
genetic algorithm to achieve state-of-the-art performance for inverse design,
with the [code here](https://github.com/aspuru-guzik-group/GA).
* SELFIES allows for [highly efficient exploration and interpolation of the chemical space](https://chemrxiv.org/articles/preprint/Beyond_Generative_Models_Superfast_Traversal_Optimization_Novelty_Exploration_and_Discovery_STONED_Algorithm_for_Molecules_using_SELFIES/13383266), with a [deterministic algorithms, see code](https://github.com/aspuru-guzik-group/stoned-selfies).
* We use SELFIES for [Deep Molecular dreaming](https://arxiv.org/abs/2012.09712), a new generative model inspired by interpretable neural networks in computational vision. See the [code of PASITHEA here](https://github.com/aspuru-guzik-group/Pasithea).
* Kohulan Rajan, Achim Zielesny, Christoph Steinbeck show in two papers that SELFIES outperforms other representations in [img2string](https://link.springer.com/article/10.1186/s13321-020-00469-w) and [string2string](https://chemrxiv.org/articles/preprint/STOUT_SMILES_to_IUPAC_Names_Using_Neural_Machine_Translation/13469202/1) translation tasks, see the codes of [DECIMER](https://github.com/Kohulan/DECIMER-Image-to-SMILES) and [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator).
* Nathan Frey, Vijay Gadepally, and Bharath Ramsundar used SELFIES with normalizing flows to develop the [FastFlows](https://arxiv.org/abs/2201.12419) framework for deep chemical generative modeling.
* An improvement to the old genetic algorithm, the authors have also released [JANUS](https://arxiv.org/abs/2106.04011), which allows for more efficient optimization in the chemical space. JANUS makes use of [STONED-SELFIES](https://pubs.rsc.org/en/content/articlepdf/2021/sc/d1sc00231g) and a neural network for efficient sampling.

## Tests
`selfies` uses `pytest` with `tox` as its testing framework.
All tests can be found in  the `tests/` directory. To run the test suite for
SELFIES, install ``tox`` and run:

```bash
tox -- --trials=10000 --dataset_samples=10000
```

By default, `selfies` is tested against a random subset
(of size ``dataset_samples=10000``) on various datasets:

 * 130K molecules from [QM9](https://www.nature.com/articles/sdata201422)
 * 250K molecules from [ZINC](https://en.wikipedia.org/wiki/ZINC_database)
 * 50K molecules from a dataset of [non-fullerene acceptors for organic solar cells](https://www.sciencedirect.com/science/article/pii/S2542435117301307)
 * 160K+ molecules from various [MoleculeNet](http://moleculenet.ai/datasets-1) datasets
 * 36M+ molecules from the [eMolecules Database](https://www.emolecules.com/info/products-data-downloads.html).
   Due to its large size, this dataset is not included on the repository. To run tests
   on it, please download the dataset into the ``tests/test_sets`` directory
   and run the ``tests/run_on_large_dataset.py`` script.

## Version History
See [CHANGELOG](https://github.com/aspuru-guzik-group/selfies/blob/master/CHANGELOG.md).

## Credits

We thank Jacques Boitreaud, Andrew Brereton, Nessa Carson (supersciencegrl), Matthew Carbone (x94carbone),  Vladimir Chupakhin (chupvl), Nathan Frey (ncfrey), Theophile Gaudin,
HelloJocelynLu, Hyunmin Kim (hmkim), Minjie Li, Vincent Mallet, Alexander Minidis (DocMinus), Kohulan Rajan (Kohulan),
Kevin Ryan (LeanAndMean), Benjamin Sanchez-Lengeling, Andrew White, Zhenpeng Yao and Adamo Young for their suggestions and bug reports,
and Robert Pollice for chemistry advices.

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
