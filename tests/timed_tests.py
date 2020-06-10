import random
import time

import selfies as sf
import selfiesv1 as sfv1


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

    print(f"Timing {len(smiles)} SMILES from {file_path}")

    # time sfv1.encoder
    start = time.time()
    for s in smiles:
        sfv1.encoder(s)
    enc_time = time.time() - start
    print(f"--> selfies.encoder: {enc_time:0.7f}s")

    # time sfv1.decoder
    start = time.time()
    for s in selfies:
        sfv1.decoder(s)
    dec_time = time.time() - start
    print(f"--> selfies.decoder: {dec_time:0.7f}s")


if __name__ == '__main__':

    # temporary example
    time_roundtrip('test_sets/dataA_QM9.txt', -1)
