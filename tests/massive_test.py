import time

import pandas as pd
from rdkit.Chem import MolFromSmiles, MolToSmiles

import selfies as sf

error_list = []
with open("test_sets/error_list.csv", "w") as error_log:
    error_log.write("In, Out, Canonical In, Canonical Out\n")

print("Starting Test")

total_decoded = 0
total_time = 0

test_path = 'test_sets/dataZ_massive.smi'

for chunk in pd.read_csv(test_path, chunksize=10000, delimiter=' '):

    for smiles in chunk['isosmiles']:

        try:

            start = time.time()
            encoded = sf.encoder(smiles)
            decoded_smiles = sf.decoder(encoded)
            total_time += (time.time() - start)

            if MolFromSmiles(decoded_smiles) is None:
                raise ValueError

            can_input = MolToSmiles(MolFromSmiles(smiles))
            can_output = MolToSmiles(MolFromSmiles(decoded_smiles))

            if can_input != can_output:
                error_list.append((smiles, decoded_smiles,
                                   can_input, can_output))

        except (Exception, ValueError):
            error_list.append((smiles, "", "", ""))

    total_decoded += 10000

    print(f"{total_decoded: 9} Molecules Decoded\t"
          f"{total_time:20.5f}s\t"
          f"{len(error_list): 9} Errors")

    with open("test_sets/error_list.csv", "a") as error_log:
        for error in error_list:
            error_log.write(','.join(error) + "\n")
    error_list = []
