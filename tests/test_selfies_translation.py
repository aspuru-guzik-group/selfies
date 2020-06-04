import faulthandler

faulthandler.enable()

import pytest
import random

import selfies as sf
from rdkit.Chem import MolFromSmiles, MolToSmiles

test_path_list = [
    'test_sets/dataA_QM9.txt',
    'test_sets/dataB_NonFullerene.txt',
    'test_sets/dataJ_250k_rndm_zinc_drugs_clean.txt'
]


@pytest.mark.parametrize("test_paths", [test_path_list])
def test_roundtrip_translation(test_paths, sample_size=1000):
    """Tests <sample_size> random SMILES from each data set in <test_paths>;
    encodes and then decodes the SMILES, and checks whether the decoded
    SMILES corresponds to the same molecule as the input SMILES.
    """

    error_list = []

    for test_path in test_paths:

        with open(test_path, 'r') as test_file:
            smiles_set = random.sample(test_file.readlines(), sample_size)
            smiles_set = [line.rstrip() for line in smiles_set]

        for smiles in smiles_set:

            decoded_smiles = sf.decoder(sf.encoder(smiles), N_restrict=False)

            can_input = MolToSmiles(MolFromSmiles(smiles))
            can_output = MolToSmiles(MolFromSmiles(decoded_smiles))

            if can_input != can_output:
                error_list.append((smiles, decoded_smiles,
                                   can_input, can_output))

    with open("error_list.csv", "w") as error_log:
        error_log.write("In, Out, Canonical In, Canonical Out")
        for error in error_list:
            error_log.write(','.join(error) + "\n")

    assert len(error_list) == 0


def test_random_selfies_decoder():
    """Passes random strings built from the SELFIES alphabet to the decoder,
    and uses RDKit to check whether the decoded SMILES are valid.
    """
    trials = 1000
    max_len = 50
    alphabet = sf.selfies_alphabet()

    for i in range(trials):

        # create random SELFIES and decode
        rand_len = random.randint(1, max_len)
        selfies = ''.join(random.choice(alphabet) for _ in range(rand_len))
        smiles = sf.decoder(selfies)

        # check if SMILES is valid
        try:
            is_valid = MolFromSmiles(smiles, sanitize=True) is not None
        except Exception:
            is_valid = False

        assert is_valid, f"Invalid SMILES {smiles} decoded from {selfies}"


if __name__ == '__main__':
    pytest.main()
