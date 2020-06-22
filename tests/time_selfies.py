import random
import time

import os
import selfies as sf

from rdkit.Chem import MolToSmiles, MolFromSmiles, Kekulize


def time_roundtrip(file_path: str, sample_size: int = -1):
    """Tests the amount of time it takes to encode and then decode an
    entire .txt file of SMILES strings <n> times. If <sample_size> is positive,
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

        for i in range(len(smiles)):

            if 'c' in smiles[i]:
                m = MolFromSmiles(smiles[i])
                Kekulize(m)
                smiles[i] = MolToSmiles(m, kekuleSmiles=True)

    print(f"Timing {len(smiles)} SMILES from {file_path}")

    # time sf.encoder
    start = time.time()
    for s in smiles:
        sf.encoder(s)
    enc_time = time.time() - start
    print(f"--> selfies.encoder: {enc_time:0.7f}s")

    # time sf.decoder
    start = time.time()
    for s in selfies:
        sf.decoder(s)
    dec_time = time.time() - start
    print(f"--> selfies.decoder: {dec_time:0.7f}s")


if __name__ == '__main__':

    # temporary example
    time_roundtrip('test_sets/130K_QM9.txt')
