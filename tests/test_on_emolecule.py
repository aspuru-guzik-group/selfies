""" This script is a clone of the 'test_on_datasets.py' script,
specifically for testing the large 1.2 GB eKolecules set. See
tests/README.md for instructions on how to download, process, and
test the set for developers interested in making a PR to the repository!

This file is split into 23 text files of 1 million SMILES strings,
such that the developer can end the test and run from the last
file tested, by passing the specific paths into 'test_sets = {}'.
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

test_sets = [
    ('test_sets/split00.txt', 'isosmiles'),
    ('test_sets/split01.txt', 'isosmiles'),
    ('test_sets/split02.txt', 'isosmiles'),
    ('test_sets/split03.txt', 'isosmiles'),
    ('test_sets/split04.txt', 'isosmiles'),
    ('test_sets/split05.txt', 'isosmiles'),
    ('test_sets/split06.txt', 'isosmiles'),
    ('test_sets/split07.txt', 'isosmiles'),
    ('test_sets/split08.txt', 'isosmiles'),
    ('test_sets/split09.txt', 'isosmiles'),
    ('test_sets/split10.txt', 'isosmiles'),
    ('test_sets/split11.txt', 'isosmiles'),
    ('test_sets/split12.txt', 'isosmiles'),
    ('test_sets/split13.txt', 'isosmiles'),
    ('test_sets/split14.txt', 'isosmiles'),
    ('test_sets/split15.txt', 'isosmiles'),
    ('test_sets/split16.txt', 'isosmiles'),
    ('test_sets/split17.txt', 'isosmiles'),
    ('test_sets/split18.txt', 'isosmiles'),
    ('test_sets/split19.txt', 'isosmiles'),
    ('test_sets/split20.txt', 'isosmiles'),
    ('test_sets/split21.txt', 'isosmiles'),
    ('test_sets/split22.txt', 'isosmiles')
]


# @pytest.mark.parametrize("test_path, column_name", test_sets)
@pytest.mark.skip(reason="no way of currently testing this")
def test_roundtrip_translation(test_path, column_name, dataset_samples):
    """Tests a roundtrip SMILES -> SELFIES -> SMILES translation of the
    SMILES examples in QM9, NonFullerene, Zinc, etc.
    """

    constraints = sf.get_semantic_constraints()
    constraints['N'] = 6
    sf.set_semantic_constraints(constraints)

    # file I/O
    test_name = os.path.splitext(os.path.basename(test_path))[0]

    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, test_path)
    error_path = os.path.join(curr_dir, f"error_sets/errors_{test_name}.csv")

    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []
    with open(error_path, "w+") as error_log:
        error_log.write("In, Out\n")
    error_found_flag = False

    # make pandas reader
    N = sum(1 for _ in open(test_path)) - 1
    S = dataset_samples if (0 < dataset_samples <= N) else N
    skip = sorted(random.sample(range(1, N + 1), N - S))
    reader = pd.read_csv(test_path,
                         chunksize=10000,
                         header=0,
                         delimiter=' ',
                         skiprows=skip)

    # roundtrip testing
    for chunk in reader:
        for in_smiles in chunk[column_name]:

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

    assert not error_found_flag


# @pytest.mark.parametrize("test_path, column_name", test_sets)
@pytest.mark.skip(reason="no way of currently testing this")
def test_kekulize_parser(test_path, column_name, dataset_samples):
    """Tests the kekulization of SMILES, which is the first step of
    selfies.encoder().
    """

    # file I/O
    test_name = os.path.splitext(os.path.basename(test_path))[0]

    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, test_path)
    error_path = os.path.join(curr_dir,
                              f"error_sets/errors_kekulize_{test_name}.csv")

    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []
    with open(error_path, "w+") as error_log:
        error_log.write("In\n")
    error_found_flag = False

    # make pandas reader
    N = sum(1 for _ in open(test_path)) - 1
    S = dataset_samples if (0 < dataset_samples <= N) else 0
    skip = sorted(random.sample(range(1, N + 1), N - S))
    reader = pd.read_csv(test_path,
                         chunksize=10000,
                         header=0,
                         delimiter=' ',
                         skiprows=skip)

    # kekulize testing
    for chunk in reader:
        for smiles in chunk[column_name]:

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
