import faulthandler
import pathlib
import random

import pandas as pd
import pytest
from rdkit import Chem

import selfies as sf

faulthandler.enable()

TEST_SET_DIR = pathlib.Path(__file__).parent / "test_sets"
ERROR_LOG_DIR = pathlib.Path(__file__).parent / "error_logs"
ERROR_LOG_DIR.mkdir(exist_ok=True, parents=True)

datasets = list(TEST_SET_DIR.glob("**/*.csv"))


@pytest.mark.parametrize("test_path", datasets)
def test_roundtrip_translation(test_path, dataset_samples):
    """Tests SMILES -> SELFIES -> SMILES translation on various datasets.
    """

    # very relaxed constraints
    constraints = sf.get_preset_constraints("hypervalent")
    constraints.update({"P": 7, "P-1": 8, "P+1": 6, "?": 12})
    sf.set_semantic_constraints(constraints)

    error_path = ERROR_LOG_DIR / "{}.csv".format(test_path.stem)
    with open(error_path, "w+") as error_log:
        error_log.write("In, Out\n")

    error_data = []
    error_found = False

    n_lines = sum(1 for _ in open(test_path)) - 1
    n_keep = dataset_samples if (0 < dataset_samples <= n_lines) else n_lines
    skip = random.sample(range(1, n_lines + 1), n_lines - n_keep)
    reader = pd.read_csv(test_path, chunksize=10000, header=0, skiprows=skip)

    for chunk in reader:

        for in_smiles in chunk["smiles"]:
            in_smiles = in_smiles.strip()

            mol = Chem.MolFromSmiles(in_smiles, sanitize=True)
            if (mol is None) or ("*" in in_smiles):
                continue

            try:
                selfies = sf.encoder(in_smiles, strict=True)
                out_smiles = sf.decoder(selfies)
            except (sf.EncoderError, sf.DecoderError):
                error_data.append((in_smiles, ""))
                continue

            if not is_same_mol(in_smiles, out_smiles):
                error_data.append((in_smiles, out_smiles))

        with open(error_path, "a") as error_log:
            for entry in error_data:
                error_log.write(",".join(entry) + "\n")

        error_found = error_found or error_data
        error_data = []

    sf.set_semantic_constraints()  # restore constraints

    assert not error_found


def is_same_mol(smiles1, smiles2):
    try:
        can_smiles1 = Chem.CanonSmiles(smiles1)
        can_smiles2 = Chem.CanonSmiles(smiles2)
        return can_smiles1 == can_smiles2
    except Exception:
        return False
