"""Script for testing selfies against large datasets.
"""

import argparse
import pathlib

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

import selfies as sf

parser = argparse.ArgumentParser()
parser.add_argument("--data_path", type=str, default="version.smi.gz")
parser.add_argument("--col_name", type=str, default="isosmiles")
parser.add_argument("--sep", type=str, default=r"\s+")
parser.add_argument("--start_from", type=int, default=0)
args = parser.parse_args()

TEST_DIR = pathlib.Path(__file__).parent
TEST_SET_PATH = TEST_DIR / "test_sets" / args.data_path
ERROR_LOG_DIR = TEST_DIR / "error_logs"
ERROR_LOG_DIR.mkdir(exist_ok=True, parents=True)


def make_reader():
    return pd.read_csv(TEST_SET_PATH, sep=args.sep, chunksize=10000)


def roundtrip_translation():
    sf.set_semantic_constraints("hypervalent")

    n_entries = 0
    for chunk in make_reader():
        n_entries += len(chunk)
    pbar = tqdm(total=n_entries)

    reader = make_reader()
    error_log = open(ERROR_LOG_DIR / f"{TEST_SET_PATH.stem}.txt", "a+")

    curr_idx = 0
    for chunk_idx, chunk in enumerate(reader):
        for in_smiles in chunk[args.col_name]:
            pbar.update(1)
            curr_idx += 1
            if curr_idx < args.start_from:
                continue

            in_smiles = in_smiles.strip()

            mol = Chem.MolFromSmiles(in_smiles, sanitize=True)
            if (mol is None) or ("*" in in_smiles):
                continue

            try:
                selfies = sf.encoder(in_smiles, strict=True)
                out_smiles = sf.decoder(selfies)
            except (sf.EncoderError, sf.DecoderError):
                error_log.write(in_smiles + "\n")
                tqdm.write(in_smiles)
                continue

            if not is_same_mol(in_smiles, out_smiles):
                error_log.write(in_smiles + "\n")
                tqdm.write(in_smiles)

    error_log.close()


def is_same_mol(smiles1, smiles2):
    try:
        can_smiles1 = Chem.CanonSmiles(smiles1)
        can_smiles2 = Chem.CanonSmiles(smiles2)
        return can_smiles1 == can_smiles2
    except Exception:
        return False


if __name__ == "__main__":
    roundtrip_translation()
