"""This script is specifically for testing the large eMolecules dataset,
which upon downloading, should have the file name version.smi.gz.

This test will automatically create a checkpoint file, so that the test
can be paused and resumed easily. Please note that there are some SMILES
in the dataset that are expected to fail; these are SMILES which use
the same ring-bond number across a dot symbol.

Example: C1.C2.C12, COC(=O)C1CCCN1CP123OCCO1.C(CO2)O3
"""

import faulthandler
import os

import pandas as pd
import pytest
from rdkit.Chem import MolFromSmiles, MolToSmiles

import selfies as sf

faulthandler.enable()

# Test Configuration ==========================================================

# file path
curr_dir = os.path.dirname(__file__)
EMOL_PATH = os.path.join(curr_dir, 'test_sets', 'version.smi.gz')

# SMILES in eMolecules are under this column name
COL_NAME = 'isosmiles'


# Tests =======================================================================

# TODO: Comment out the pytest skip to use this script.
@pytest.mark.skip(reason="eMolecules dataset not on GitHub")
def test_roundtrip_translation():
    """Tests a roundtrip SMILES -> SELFIES -> SMILES translation of the
    SMILES examples in QM9, NonFullerene, Zinc, etc.
    """

    # modify constraints
    constraints = sf.get_semantic_constraints()
    constraints['N'] = 6
    constraints['Br'] = 7
    constraints['Cl'] = 7
    constraints['I'] = 7
    sf.set_semantic_constraints(constraints)

    # file I/O
    ckpt_path = os.path.join(curr_dir, 'checkpoints', 'emolecule_ckpt.txt')
    error_path = os.path.join(curr_dir, 'error_sets', 'errors_emolecules.csv')

    # check if a previous checkpoint exists to continue tests
    if os.path.exists(ckpt_path):
        with open(ckpt_path, 'r') as ckpt_file:
            checkpoint = int(ckpt_file.readlines()[0])

    # if no path to a checkpoint exists, create a new directory for error logging and checkpoints
    else:
        os.makedirs(os.path.dirname(ckpt_path), exist_ok=True)
        os.makedirs(os.path.dirname(error_path), exist_ok=True)

        with open(error_path, "w+") as error_log:
            error_log.write("In, Out\n")
        checkpoint = -1

    error_list = []
    error_found_flag = False

    # make pandas reader
    reader = pd.read_csv(EMOL_PATH,
                         chunksize=10000,
                         compression='gzip',
                         delimiter=' ',
                         header=0)

    # roundtrip testing
    for chunk_idx, chunk in enumerate(reader):

        if chunk_idx <= checkpoint:
            continue

        for in_smiles in chunk[COL_NAME]:

            # check if SMILES in chunk is a valid RDKit molecule. if not, skip testing
            # All inputted SMILES must be valid RDKit Mol objects to be encoded.
            if (MolFromSmiles(in_smiles) is None) or ('*' in in_smiles):
                continue

            # encode selfies
            selfies = sf.encoder(in_smiles)

            # if unable to encode SMILES, write to list of errors
            if selfies is None:
                error_list.append((in_smiles, ''))
                continue

            # take encoeded SELFIES and decode
            out_smiles = sf.decoder(selfies)

            # compare original SMILES to decoded SELFIE string, if wrong, write to list of errors.
            if not is_same_mol(in_smiles, out_smiles):
                error_list.append((in_smiles, out_smiles))

        # open and write all errors to errors_emolecule.csv
        with open(error_path, "a") as error_log:
            for error in error_list:
                error_log.write(','.join(error) + "\n")
        error_found_flag = error_found_flag or error_list
        error_list = []

        # create checkpoint from the current pandas reader chunk, to load from and continue testing.
        with open(ckpt_path, 'w+') as ckpt_file:
            ckpt_file.write(str(chunk_idx))

    sf.set_semantic_constraints()  # restore defaults
    os.remove(ckpt_path)  # remove checkpoint

    assert not error_found_flag


# Helper Methods

def is_same_mol(smiles1, smiles2):
    """Helper method that returns True if smiles1 and smiles2 correspond
    to the same molecule.
    """

    if smiles1 is None or smiles2 is None:
        return False

    m1 = MolFromSmiles(smiles1)
    m2 = MolFromSmiles(smiles2)

    if m1 is None or m2 is None:
        return False

    can1 = MolToSmiles(m1)
    can2 = MolToSmiles(m2)

    return can1 == can2
