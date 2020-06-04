import faulthandler
faulthandler.enable()

from selfies import encoder, decoder
from rdkit.Chem import MolFromSmiles, MolToSmiles
import pandas as pd
import random

test_path_list = [
    'test_sets/dataA_QM9.txt',
    'test_sets/dataB_NonFullerene.txt',
    'test_sets/dataJ_250k_rndm_zinc_drugs_clean.txt'
]

def test_selfies_translation(test_path, sample_size):
    error_list = []


    print ('Testing following dataset:' + test_path)
    smi_df = pd.read_csv(test_path)
    smi_df.columns = ['smiles']
    smiles = smi_df['smiles']
    # random.shuffle(smiles) # shuffle the SMILES such that new batch is always tested

    for i, smi in enumerate(smiles):
        if i==sample_size+1:
            break

        encoded_selfies=encoder(smi)
        decoded_smiles=decoder(encoded_selfies, N_restrict=False)
        s1=MolToSmiles(MolFromSmiles(smi))
        s2=MolToSmiles(MolFromSmiles(decoded_smiles))

        if s1==s2:
            # print('SMILES are the same: ',s1)
            continue
        else:
            print('ERROR: SMILES are different: ' + s1 + ' - ' + s2)
            error_list.append(str(i) + ': ERROR!: ' + s1)

    with open("error_list.txt", "w") as file:
        for error in error_list:
            file.write(error + "\n")

if __name__ == '__main__':
    for test_path in test_path_list:
        test_selfies_translation(test_path, 1000)
