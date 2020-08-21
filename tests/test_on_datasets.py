"""Lengthy tests that are run on testing data sets.
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

datasets = [
    ('Custom_Cases', 'smiles'),
    ('130K_QM9', 'smiles'),
    ('51K_NonFullerene', 'smiles'),
    ('250K_ZINC', 'smiles'),
    ('8k_Tox21', 'smiles'),
    ('93k_PubChem_MUV_bioassay', 'smiles')
]


@pytest.mark.parametrize("test_name, column_name", datasets)
def test_roundtrip_translation(test_name, column_name, dataset_samples):
    """Tests a roundtrip SMILES -> SELFIES -> SMILES translation of the
    SMILES examples in QM9, NonFullerene, Zinc, etc.
    """

    # modify semantic bond constraints
    constraints = sf.get_semantic_constraints()
    constraints['N'] = 6
    sf.set_semantic_constraints(constraints)

    # file I/O
    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, 'test_sets', test_name + ".txt")
    error_path = os.path.join(curr_dir,
                              'error_sets',
                              "errors_{}.csv".format(test_name))

    # create error directory
    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []

    # add header in error log text file
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
                         skiprows=skip)

    # roundtrip testing
    for chunk in reader:
        for in_smiles in chunk[column_name]:
            # check if SMILES in chunk is a valid RDKit molecule.
            # if not, skip testing
            # All inputted SMILES must be valid
            # RDKit Mol objects to be encoded.
            if (MolFromSmiles(in_smiles) is None) or ('*' in in_smiles):
                continue

            # encode SELFIE string
            selfies = sf.encoder(in_smiles)

            # if unable to encode SMILES, write to list of errors
            if selfies is None:
                error_list.append((in_smiles, ''))
                continue

            # take encoeded SELFIES and decode
            out_smiles = sf.decoder(selfies)

            # compare original SMILES to decoded SELFIE string.
            # if not the same string, write to list of errors.
            if not is_same_mol(in_smiles, out_smiles):
                error_list.append((in_smiles, out_smiles))

        # open and write all errors to errors_{test_name}.csv
        with open(error_path, "a") as error_log:
            for error in error_list:
                error_log.write(','.join(error) + "\n")
        error_found_flag = error_found_flag or error_list
        error_list = []

    sf.set_semantic_constraints()  # restore defaults

    assert not error_found_flag


@pytest.mark.skip(reason="covered by round-trip test")
@pytest.mark.parametrize("test_name, column_name", datasets)
def test_kekulize_parser(test_name, column_name, dataset_samples):
    """Tests the kekulization of SMILES, which is the first step of
    selfies.encoder().
    """

    # file I/O
    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, 'test_sets', test_name + ".txt")
    error_path = os.path.join(curr_dir,
                              'error_sets',
                              "errors_kekulize_{}.csv".format(test_name))

    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []
    with open(error_path, "w+") as error_log:
        error_log.write("In\n")
    error_found_flag = False

    # make pandas reader
    N = sum(1 for _ in open(test_path)) - 1
    S = dataset_samples if (0 < dataset_samples <= N) else N
    skip = sorted(random.sample(range(1, N + 1), N - S))
    reader = pd.read_csv(test_path,
                         chunksize=10000,
                         header=0,
                         skiprows=skip)

    # kekulize testing
    for chunk in reader:
        for smiles in chunk[column_name]:

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
