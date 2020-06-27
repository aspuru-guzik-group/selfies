import os
import time

import pandas as pd
import pytest
from rdkit.Chem import MolFromSmiles, MolToSmiles

from selfies.encoder import _parse_smiles
from selfies.kekulize import BRANCH_TYPE, RING_TYPE, kekulize_parser

test_sets = [
    ('test_sets/130K_QM9.txt', 'smiles'),
    ('test_sets/51K_NonFullerene.txt', 'smiles'),
    ('test_sets/250k_ZINC.txt', 'smiles')
]  # add if desired ('22M_eMolecule.smi', 'isosmiles')


@pytest.mark.parametrize("test_path, column_name", test_sets)
def test_kekulize_parser(test_path, column_name):
    # file I/O
    test_name = os.path.basename(test_path)
    test_name = os.path.splitext(test_name)[0]

    curr_dir = os.path.dirname(__file__)
    test_path = os.path.join(curr_dir, test_path)
    error_path = os.path.join(curr_dir,
                              f"error_sets/errors_kekulize_{test_name}.csv")

    os.makedirs(os.path.dirname(error_path), exist_ok=True)
    error_list = []
    with open(error_path, "w+") as error_log:
        error_log.write("In\n")
    error_found_flag = False

    elapsed = 0

    # roundtrip testing
    for chunk in pd.read_csv(test_path, chunksize=10000, delimiter=' '):

        for smiles in chunk[column_name]:

            try:
                start = time.time()
                kekule_gen = list(kekulize_parser(_parse_smiles(smiles)))
                elapsed += time.time() - start

                kekule_smiles = ""
                for bond, char, char_type in kekule_gen:
                    if char_type == BRANCH_TYPE:
                        bond = ''
                    kekule_smiles += bond

                    if char_type == RING_TYPE and len(char) == 2:
                        kekule_smiles += "%"
                    kekule_smiles += char

                can_smiles = MolToSmiles(MolFromSmiles(smiles))
                can_kekule = MolToSmiles(MolFromSmiles(kekule_smiles,
                                                       sanitize=True))

                if can_smiles != can_kekule:
                    raise ValueError

            except (ValueError, Exception) as err:
                print(err, smiles)
                error_list.append(smiles)

        with open(error_path, "a") as error_log:
            error_log.write("\n".join(error_list))
        error_found_flag = error_found_flag or error_list
        error_list = []

    print()
    print(f"{test_name} in {elapsed:11.7f}s")
    assert not error_found_flag
