"""This script is specifically for testing the large eMolecules dataset.
See tests/README.md for instructions on how to download, process, and
test the set for developers interested in making a PR to the repository!

This file is split into text files of 1 million SMILES strings,
such that the developer can end the test and run from the last
file tested, by passing the specific paths into 'datasets = []'.
"""

import faulthandler
import os
import random

import pandas as pd
import pytest
from rdkit.Chem import MolFromSmiles, MolToSmiles

import selfies as sf
from selfies.encoder import _parse_smiles
from selfies.kekulize import BRANCH_TYPE, RING_TYPE, kekulize_parser

faulthandler.enable()

# Test Configuration ==========================================================
# TODO: Edit the configurations of the eMolecule script below.

datasets = []

# SMILES in eMolecules are under this column name
COL_NAME = 'isosmiles'

# How many samples from each data file to test. -1 to test the whole file.
NUM_SAMPLES = -1

# Dynamically search for all text files in the eMolecules folder
curr_dir = os.path.dirname(__file__)
emol_dir = os.path.join(curr_dir, 'test_sets', 'eMolecules')

for file in os.listdir(emol_dir):
    if file.endswith(".txt"):
        datasets.append(file[:-4])  # remove the .txt extension


# Alternatively, specify the files you want to test:
# e.g. datasets = ['split05', 'split06']


# Tests =======================================================================

# TODO: Comment out the pytest skip to use this script.
@pytest.mark.skip(reason="eMolecules dataset not on GitHub")
@pytest.mark.parametrize("test_name", datasets)
def test_roundtrip_translation(test_name):
    """Tests a roundtrip SMILES -> SELFIES -> SMILES translation of the
    SMILES examples in QM9, NonFullerene, Zinc, etc.
    """

    constraints = sf.get_semantic_constraints()
    constraints['N'] = 6
    constraints['Br'] = 7
    constraints['Cl'] = 7
    constraints['I'] = 7
    sf.set_semantic_constraints(constraints)

    # file I/O
    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, 'test_sets', 'eMolecules',
                             f"{test_name}.txt")
    error_path = os.path.join(curr_dir,
                              'error_sets',
                              f"errors_{test_name}.csv")

    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []
    with open(error_path, "w+") as error_log:
        error_log.write("In, Out\n")
    error_found_flag = False

    # make pandas reader
    N = sum(1 for _ in open(test_path)) - 1
    S = NUM_SAMPLES if (0 < NUM_SAMPLES <= N) else N
    skip = sorted(random.sample(range(1, N + 1), N - S))
    reader = pd.read_csv(test_path,
                         chunksize=10000,
                         header=0,
                         skiprows=skip)

    # roundtrip testing
    for chunk in reader:
        for in_smiles in chunk[COL_NAME]:

            if (MolFromSmiles(in_smiles) is None) or ('*' in in_smiles):
                continue

            selfies = sf.encoder(in_smiles)

            if selfies is None:
                error_list.append((in_smiles, ''))
                continue

            out_smiles = sf.decoder(selfies)

            if not is_same_mol(in_smiles, out_smiles):
                error_list.append((in_smiles, out_smiles))

        with open(error_path, "a") as error_log:
            for error in error_list:
                error_log.write(','.join(error) + "\n")
        error_found_flag = error_found_flag or error_list
        error_list = []

    sf.set_semantic_constraints()  # restore defaults

    assert not error_found_flag


# TODO: Comment out the pytest skip to use this script.
#  This test is somewhat covered by the first test though.
@pytest.mark.skip(reason="eMolecules dataset not on GitHub")
@pytest.mark.parametrize("test_name", datasets)
def test_kekulize_parser(test_name):
    """Tests the kekulization of SMILES, which is the first step of
    selfies.encoder().
    """

    # file I/O
    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, 'test_sets', 'eMolecules',
                             f"{test_name}.txt")
    error_path = os.path.join(curr_dir,
                              'error_sets',
                              f"errors_kekulize_{test_name}.csv")

    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []
    with open(error_path, "w+") as error_log:
        error_log.write("In\n")
    error_found_flag = False

    # make pandas reader
    N = sum(1 for _ in open(test_path)) - 1
    S = NUM_SAMPLES if (0 < NUM_SAMPLES <= N) else N
    skip = sorted(random.sample(range(1, N + 1), N - S))
    reader = pd.read_csv(test_path,
                         chunksize=10000,
                         header=0,
                         skiprows=skip)

    # kekulize testing
    for chunk in reader:
        for smiles in chunk[COL_NAME]:

            if (MolFromSmiles(smiles) is None) or ('*' in smiles):
                continue

            # build kekulized SMILES
            kekule_fragments = []

            for fragment in smiles.split("."):

                kekule_gen = kekulize_parser(_parse_smiles(fragment))

                k = []
                for bond, symbol, symbol_type in kekule_gen:
                    if symbol_type == BRANCH_TYPE:
                        bond = ''
                    k.append(bond)

                    if symbol_type == RING_TYPE and len(symbol) == 2:
                        k.append('%')
                    k.append(symbol)

                kekule_fragments.append(''.join(k))

            kekule_smiles = '.'.join(kekule_fragments)

            if not is_same_mol(smiles, kekule_smiles):
                error_list.append(smiles)

        with open(error_path, "a") as error_log:
            error_log.write("\n".join(error_list))
        error_found_flag = error_found_flag or error_list
        error_list = []

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
