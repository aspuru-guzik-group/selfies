This directory contains the testing suite for ``selfies``. 
 * ``test_sets/`` contains the testing datasets.
 * ``test_on_datasets.py`` contains tests that use the testing datasets.
 * ``test_on_emolecules.py`` contains tests on the eMolecules Plus dataset. 
    Due to its large size, this dataset is not included on the repository. To run tests 
    on it, please download the dataset in the ``tests/test_sets`` directory 
    and enable its pytest at ``tests/test_on_emolecules.py``. 
 * ``test_selfies.py`` contains tests for the overall ``selfies`` library. 
 * ``test_selfies_utils`` contains tests for the utility methods 
    of ``selfies``, for example ``selfies.len_selfies``.
 * ``time_selfies.py`` contains a script that times the efficiency of 
    ``selfies``.
