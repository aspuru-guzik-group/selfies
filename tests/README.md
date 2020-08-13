This directory contains the testing suite for ``selfies``. 
 * ``test_sets/`` contains the testing datasets.
 * ``test_on_datasets.py`` contains tests that use the testing datasets.
 * ``test_on_emolecule.py`` contains tests on the eMolecules Plus dataset. 
    The dataset must be downloaded and processed, for which instructions are 
    included below and in the main README.
 * ``test_selfies.py`` contains tests for the overall ``selfies`` library. 
 * ``test_selfies_utils`` contains tests for the utility methods 
    of ``selfies``, for example ``selfies.len_selfies``.
 * ``time_selfies.py`` contains a script that times the efficiency of 
    ``selfies``.

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

