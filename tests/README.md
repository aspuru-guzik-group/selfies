This directory contains the testing suite for ``selfies``. 
 * ``test_sets/`` contains the testing datasets.
 * ``test_on_datasets.py`` contains tests that use the testing datasets.
 * ``test_on_emolecule.py`` contains tests on the eMolecules dataset of over 22 million SMILES strings. The dataset must be downloaded and processed, for which instructions are included below.
 * ``test_selfies.py`` contains tests for the overall ``selfies`` library. 
 * ``test_selfies_utils`` contains tests for the utility methods 
    of ``selfies``, for example ``selfies.len_selfies``.
 * ``time_selfies.py`` contains a script that times the efficiency of 
    ``selfies``.

---

The eMolecules dataset is not included in the github repository, due to its size.
For developers interested on testing using this set, one can download the 
dataset [here](https://www.emolecules.com/info/plus/download-database). 
Once the dataset is downloaded, one can split the file into chunks of 1 
million SMILES for easier processing using the following command in Windows:

```
split -l 1000000 -d --additional-suffix=.txt version.txt split
```

On Mac, the following commands can be used:
```
brew install coreutils
gsplit -l 1000000 -d --additional-suffix=.txt version.txt split
```

Once the dataset is downloaded and split in the `test_sets/` directory,
 `test_on_emolecule.py` can be run to test the eMolecules set. 

