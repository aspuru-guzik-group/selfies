from rdkit import Chem

import selfies as sf


def test_oversized_ring():
    """Test SELFIES that have a ring, with Q so large that the (Q + 1)-th
    previously derived atom does not exist.
    """

    assert is_eq(sf.decoder("[C][C][C][C][Ring1][O]"), "C1CCC1")


def test_ring_at_beginning_of_branch():
    """Test SELFIES that have a ring immediately at the start of a branch.
    """

    assert is_eq(sf.decoder("[C][C][C][C][Branch1_1][Branch1_1]"
                            "[Ring1][Ring2][C][Cl][F]"),
                 "C1CCC1(CCl)F")


def test_ring_imediately_following_branch():
    """Test SELFIES that have a ring immediately following after a branch.
    """

    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[Ring1][Branch1_1]"),
                 "CCC1CCCC1(OCO)")
    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[Branch1_1][C][F][Ring1][Branch1_1]"),
                 "CCC1CCCC1(OCO)(F)")


def test_ring_after_branch():
    """Tests SELFIES that have a ring following a branch, but not
    immediately after a branch.
    """

    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[C][Ring1][Branch1_1]"),
                 "CCCCCCC(OCO)=C")
    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[Branch1_1][C][F][C][C][Ring1][Branch2_2]"),
                 "CCCCC1CC(OCO)(F)CC1")


def test_ring_at_end_of_selfies():
    """Test SELFIES that have a ring symbol as its very last character.
    """

    assert is_eq(sf.decoder("[C][C][C][C][C][Ring1]"), "CCCC=C")


def test_ring_on_top_of_existing_bond():
    """Tests SELFIES with rings between two atoms that are already bonded
    in the main scaffold.
    """

    assert is_eq(sf.decoder("[C][C][Ring1][C]"), "C=C")
    assert is_eq(sf.decoder("[C][/C][Ring1][C]"), "C=C")
    assert is_eq(sf.decoder("[C][C][Expl=Ring1][C]"), "C#C")
    assert is_eq(sf.decoder("[C][C][Expl#Ring1][C]"), "C#C")


def test_consecutive_rings():
    """Test SELFIES which have multiple consecutive rings.
    """

    assert is_eq(sf.decoder("[C][C][C][C][Ring1][Ring2][Ring1][Ring2]"),
                 "C=1CCC=1")  # 1 + 1
    assert is_eq(sf.decoder("[C][C][C][C][Ring1][Ring2][Ring1][Ring2]"
                            "[Ring1][Ring2]"),
                 "C#1CCC#1")  # 1 + 1 + 1
    assert is_eq(sf.decoder("[C][C][C][C][Expl=Ring1][Ring2][Ring1][Ring2]"),
                 "C#1CCC#1")  # 2 + 1
    assert is_eq(sf.decoder("[C][C][C][C][Ring1][Ring2][Expl=Ring1][Ring2]"),
                 "C#1CCC#1")  # 1 + 2

    # consecutive rings that exceed bond constraints
    assert is_eq(sf.decoder("[C][C][C][C][Expl#Ring1][Ring2]"
                            "[Expl=Ring1][Ring2]"),
                 "C#1CCC#1")  # 3 + 2
    assert is_eq(sf.decoder("[C][C][C][C][Expl=Ring1][Ring2]"
                            "[Expl#Ring1][Ring2]"),
                 "C#1CCC#1")  # 2 + 3
    assert is_eq(sf.decoder("[C][C][C][C][Expl=Ring1][Ring2]"
                            "[Expl=Ring1][Ring2]"),
                 "C#1CCC#1")  # 2 + 2

    # consecutive rings with stereochemical single bond
    assert sf.decoder("[C][C][C][C][Expl/Ring1][Ring2]") == "C/1CCC/1"
    assert sf.decoder("[C][C][C][C][Expl/Ring1][Ring2][Ring1][Ring2]") \
           == "C=1CCC=1"


def test_unconstrained_symbol():
    """Tests SELFIES with symbols that are not semantically constrained.
    """

    assert sf.decoder("[Xe-2expl][Branch1_1][C][F][Branch1_1][C][F]"
                      "[Branch1_1][C][F][Branch1_1][C][F][Branch1_1][C][F]"
                      "[Branch1_1][C][F][Branch1_1][C][F][Branch1_1][C][F]") \
           == "[Xe-2](F)(F)(F)(F)(F)(F)(F)CF"

    # change default semantic constraints
    constraints = sf.get_semantic_constraints()
    constraints['?'] = 2
    sf.set_semantic_constraints(constraints)

    assert sf.decoder("[Xe-2expl][Branch1_1][C][F][Branch1_1][C][F]"
                      "[Branch1_1][C][F][Branch1_1][C][F][Branch1_1][C][F]"
                      "[Branch1_1][C][F][Branch1_1][C][F][Branch1_1][C][F]") \
           == "[Xe-2](F)CF"

    sf.set_semantic_constraints()


# Helper Methods

def is_eq(smiles1, smiles2):
    return Chem.CanonSmiles(smiles1) == Chem.CanonSmiles(smiles2)
