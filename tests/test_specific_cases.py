from rdkit import Chem

import selfies as sf


def test_branch_and_ring_at_state_X0():
    """Tests SELFIES with branches and rings at state X0 (i.e. at the
    very beginning of a SELFIES). These symbols should be skipped.
    """

    assert is_eq(sf.decoder("[Branch3_1][C][S][C][O]"), "CSCO")
    assert is_eq(sf.decoder("[Ring3][C][S][C][O]"), "CSCO")
    assert is_eq(sf.decoder("[Branch1_1][Ring1][Ring3][C][S][C][O]"), "CSCO")


def test_branch_at_state_X1():
    """Test SELFIES with branches at state X1 (i.e. at an atom that
    can only make one bond. In this case, the branch symbol should be skipped.
    """
    assert is_eq(sf.decoder("[C][C][O][Branch1_1][C][I]"), "CCOCI")
    assert is_eq(sf.decoder("[C][C][C][O][Branch3_3][C][I]"), "CCCOCI")


def test_branch_at_end_of_selfies():
    """Test SELFIES that have a branch symbol as its very last symbol.
    """

    assert is_eq(sf.decoder("[C][C][C][C][Branch1_1]"), "CCCC")
    assert is_eq(sf.decoder("[C][C][C][C][Branch3_3]"), "CCCC")


def test_branch_with_no_atoms():
    """Test SELFIES that have a branch, but the branch has no atoms in it.
    Such branches should not be made.
    """

    assert is_eq(sf.decoder("[C][Branch1_1][Ring2][Branch1_1]"
                            "[Branch1_1][Branch1_1][F]"),
                 "CF")
    assert is_eq(sf.decoder("[C][Branch1_1][Ring2][Ring1]"
                            "[Ring1][Branch1_1][F]"),
                 "CF")

    # special case: Branch3_3 takes Q_1, Q_2 = [O] and Q_3 = ''. However,
    # there are no more symbols in the branch.
    assert is_eq(sf.decoder("[C][C][C][C][Branch3_3][O][O]"), "CCCC")


def test_oversized_ring():
    """Test SELFIES that have a ring, with Q so large that the (Q + 1)-th
    previously derived atom does not exist.
    """

    assert is_eq(sf.decoder("[C][C][C][C][Ring1][O]"), "C1CCC1")
    assert is_eq(sf.decoder("[C][C][C][C][Ring2][O][C]"), "C1CCC1")

    # special case: Ring2 takes Q_1 = [O] and Q_2 = '', leading to
    # Q = 9 * 16 + 0 (i.e. an oversized ring)
    assert is_eq(sf.decoder("[C][C][C][C][Ring2][O]"), "C1CCC1")

    # special case: ring between 1st atom and 1st atom should not be formed
    assert is_eq(sf.decoder("[C][Ring1][O]"), "C")


def test_ring_at_beginning_of_branch():
    """Test SELFIES that have a ring immediately at the start of a branch.
    """

    assert is_eq(sf.decoder("[C][C][C][C][Branch1_1][Branch1_1]"
                            "[Ring1][Ring2][C][Cl][F]"),
                 "C1CCC1(CCl)F")
    assert is_eq(sf.decoder("[C][C][C][S][Branch1_1][C][Br]"
                            "[Branch1_1][Branch1_1][Ring1][Ring2][C][Cl][F]"),
                 "C1CCS1(Br)(CCl)F")


def test_ring_immediately_following_branch():
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
    """Test SELFIES that have a ring symbol as its very last symbol.
    """

    assert is_eq(sf.decoder("[C][C][C][C][C][Ring1]"), "CCCC=C")
    assert is_eq(sf.decoder("[C][C][C][C][C][Ring3]"), "CCCC=C")


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
