import random
import time

import selfies as sf

from rdkit.Chem import MolToSmiles, MolFromSmiles, Kekulize


def time_roundtrip(file_path: str, sample_size: int = -1):
    """Tests the amount of time it takes to encode and then decode an
    entire .txt file of SMILES strings <n> times. If <sample_size> is positive,
    then a random sample is taken from the file instead.
    """

    # load data
    with open(file_path, 'r') as file:
        smiles = [line.rstrip() for line in file.readlines()]

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
    time_roundtrip('test_sets/dataA_QM9.txt', -1)
    # time_roundtrip('test_sets/dataB_NonFullerene.txt', -1)
    # time_roundtrip('test_sets/dataJ_250k_rndm_zinc_drugs_clean.txt', -1)

