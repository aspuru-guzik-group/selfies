import selfies as sf


def test_encoder_attribution():
    smiles = "C1([O-])C=CC=C1Cl"
    indices = [0, 3, 3, 3, 5, 7, 8, 10, None, None, 12]
    s, am = sf.encoder(smiles, attribute=True)
    # check that Cl lined up
    for i, ta in enumerate(am):
        if ta[1]:
            assert indices[i] == ta[1][0][0], f'found {ta[1]}; should be {indices[i]}'
        if ta[0] == '[Cl]':
            for i, v in ta[1]:
                if v == 'Cl':
                    return
    raise ValueError('Failed to find Cl in attribution map')


test_encoder_attribution()
