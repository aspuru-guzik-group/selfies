import os
import random
import time

from rdkit.Chem import Kekulize, MolFromSmiles, MolToSmiles

import selfies as sf
from selfies.encoder import _parse_smiles
from selfies.kekulize import kekulize_parser


def time_roundtrip(file_path: str, sample_size: int = -1):
    """Tests the amount of time it takes to encode and then decode an
    entire .txt file of SMILES strings. If <sample_size> is positive,
    then a random sample is taken from the file instead.
    """

    curr_dir = os.path.dirname(__file__)
    file_path = os.path.join(curr_dir, file_path)

    # load data
    with open(file_path, 'r') as file:
        smiles = [line.rstrip() for line in file.readlines()]
        smiles.pop(0)

        if sample_size > 0:
            smiles = random.sample(smiles, sample_size)
        selfies = list(map(sf.encoder, smiles))

    print("Timing {} SMILES from {}".format(len(smiles), file_path))

    # time sf.encoder
    start = time.time()
    for s in smiles:
        sf.encoder(s)
    enc_time = time.time() - start
    print("--> selfies.encoder: {:0.7f}s".format(enc_time))

    # time sf.decoder
    start = time.time()
    for s in selfies:
        sf.decoder(s)
    dec_time = time.time() - start
    print("--> selfies.decoder: {:0.7f}s".format(dec_time))


def time_kekulize(file_path: str, sample_size: int = -1):
    curr_dir = os.path.dirname(__file__)
    file_path = os.path.join(curr_dir, file_path)

    # load data
    with open(file_path, 'r') as file:
        smiles = [line.rstrip() for line in file.readlines()]
        smiles.pop(0)

        if sample_size > 0:
            smiles = random.sample(smiles, sample_size)

    print("Timing Kekulization of {} SMILES from {}".format(len(smiles),
                                                            file_path))

    # time selfies kekulization
    start = time.time()
    for s in smiles:
        list(kekulize_parser(_parse_smiles(s)))
    selfies_time = time.time() - start
    print("--> selfies kekulize: {:0.7f}s".format(selfies_time))

    # time RDKit kekulization
    start = time.time()
    for s in smiles:
        m = MolFromSmiles(s)
        Kekulize(m)
        MolToSmiles(m, kekuleSmiles=True)
    rdkit_time = time.time() - start
    print("--> RDKit kekulize: {:0.7f}s".format(rdkit_time))


if __name__ == '__main__':

    # temporary example
    time_roundtrip('test_sets/250K_ZINC.txt')
