from rdkit import Chem

import selfies as sf


def test_ring_at_beginning_of_bond():
    assert is_eq_dec("[C][C][C][C][Branch1_1][Branch1_1]"
                     "[Ring1][Ring2][C][Cl][F]",
                     "C1CCC1(CCl)F")


def test_ring_on_top_of_existing_bond():
    assert is_eq_dec("[C][C][Ring1][C]", "C=C")
    assert is_eq_dec("[C][C][Expl=Ring1][C]", "C#C")
    assert is_eq_dec("[C][C][Expl#Ring1][C]", "C#C")


def test_consecutive_rings():
    assert is_eq_dec("[C][C][C][C][Ring1][Ring2][Ring1][Ring2]",
                     "C=1CCC=1")  # 1 + 1
    assert is_eq_dec("[C][C][C][C][Ring1][Ring2][Ring1][Ring2][Ring1][Ring2]",
                     "C#1CCC#1")  # 1 + 1 + 1
    assert is_eq_dec("[C][C][C][C][Expl=Ring1][Ring2][Ring1][Ring2]",
                     "C#1CCC#1")  # 2 + 1
    assert is_eq_dec("[C][C][C][C][Ring1][Ring2][Expl=Ring1][Ring2]",
                     "C#1CCC#1")  # 1 + 2

    assert is_eq_dec("[C][C][C][C][Expl#Ring1][Ring2][Expl=Ring1][Ring2]",
                     "C#1CCC#1")  # 3 + 2
    assert is_eq_dec("[C][C][C][C][Expl=Ring1][Ring2][Expl#Ring1][Ring2]",
                     "C#1CCC#1")  # 2 + 3
    assert is_eq_dec("[C][C][C][C][Expl=Ring1][Ring2][Expl=Ring1][Ring2]",
                     "C#1CCC#1")  # 2 + 2


def is_eq_dec(selfies, smiles):
    return Chem.CanonSmiles(sf.decoder(selfies)) == Chem.CanonSmiles(smiles)
