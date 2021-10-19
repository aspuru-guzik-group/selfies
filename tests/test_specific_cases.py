import pytest

import selfies as sf


def decode_eq(selfies, smiles):
    s = sf.decoder(selfies)
    return s == smiles


def test_branch_and_ring_at_state_X0():
    """Tests SELFIES with branches and rings at state X0 (i.e. at the
    very beginning of a SELFIES). These symbols should be skipped.
    """

    assert decode_eq("[Branch3][C][S][C][O]", "CSCO")
    assert decode_eq("[Ring3][C][S][C][O]", "CSCO")
    assert decode_eq("[Branch1][Ring1][Ring3][C][S][C][O]", "CSCO")


def test_branch_at_state_X1():
    """Test SELFIES with branches at state X1 (i.e. at an atom that
    can only make one bond. In this case, the branch symbol should be skipped.
    """

    assert decode_eq("[C][C][O][Branch1][C][I]", "CCOCI")
    assert decode_eq("[C][C][C][O][#Branch3][C][I]", "CCCOCI")


def test_branch_and_ring_decrement_state():
    """Tests that the branch and ring symbols properly decrement the
    derivation state.
    """

    assert decode_eq("[C][C][C][Ring1][Ring1][#C]", "C1CC1=C")
    assert decode_eq("[C][=C][C][C][#Ring1][Ring1][#C]", "C=C1CC1")
    assert decode_eq("[C][O][C][C][=Ring1][Ring1][#C]", "COCCC")

    assert decode_eq("[C][=C][Branch1][C][=C][#C]", "C=C(C)C")


def test_branch_at_end_of_selfies():
    """Test SELFIES that have a branch symbol as its very last symbol.
    """

    assert decode_eq("[C][C][C][C][Branch1]", "CCCC")
    assert decode_eq("[C][C][C][C][#Branch3]", "CCCC")


def test_ring_at_end_of_selfies():
    """Test SELFIES that have a ring symbol as its very last symbol.
    """

    assert decode_eq("[C][C][C][C][C][Ring1]", "CCCC=C")
    assert decode_eq("[C][C][C][C][C][Ring3]", "CCCC=C")


def test_branch_with_no_atoms():
    """Test SELFIES that have a branch, but the branch has no atoms in it.
    Such branches should not be made in the outputted SMILES.
    """

    s = "[C][Branch1][Ring2][Branch1][Branch1][Branch1][F]"
    assert decode_eq(s, "CF")

    s = "[C][Branch1][Ring2][Ring1][Ring1][Branch1][F]"
    assert decode_eq(s, "CF")

    s = "[C][=Branch1][Ring2][Branch1][C][Cl][F]"
    assert decode_eq(s, "C(Cl)F")

    # special case: #Branch3 takes Q_1, Q_2 = [O] and Q_3 = ''. However,
    # there are no more symbols in the branch.
    assert decode_eq("[C][C][C][C][#Branch3][O][O]", "CCCC")


def test_oversized_branch():
    """Test SELFIES that have a branch, with Q larger than the length
    of the SELFIES
    """

    assert decode_eq("[C][Branch2][O][O][C][C][S][F][C]", "CCCSF")
    assert decode_eq("[C][#Branch2][O][O][#C][C][S][F]", "C#CCSF")


def test_oversized_ring():
    """Test SELFIES that have a ring, with Q so large that the (Q + 1)-th
    previously derived atom does not exist.
    """

    assert decode_eq("[C][C][C][C][Ring1][O]", "C1CCC1")
    assert decode_eq("[C][C][C][C][Ring2][O][C]", "C1CCC1")

    # special case: Ring2 takes Q_1 = [O] and Q_2 = '', leading to
    # Q = 9 * 16 + 0 (i.e. an oversized ring)
    assert decode_eq("[C][C][C][C][Ring2][O]", "C1CCC1")

    # special case: ring between 1st atom and 1st atom should not be formed
    assert decode_eq("[C][Ring1][O]", "C")


def test_branch_at_beginning_of_branch():
    """Test SELFIES that have a branch immediately at the start of a branch.
    """

    # [C@]((Br)Cl)F
    s = "[C@][=Branch1][Branch1][Branch1][C][Br][Cl][F]"
    assert decode_eq(s, "[C@](Br)(Cl)F")

    # [C@](((Br)Cl)I)F
    s = "[C@][#Branch1][Branch2][=Branch1][Branch1][Branch1][C][Br][Cl][I][F]"
    assert decode_eq(s, "[C@](Br)(Cl)(I)F")

    # [C@]((Br)(Cl)I)F
    s = "[C@][#Branch1][Branch2][Branch1][C][Br][Branch1][C][Cl][I][F]"
    assert decode_eq(s, "[C@](Br)(Cl)(I)F")


def test_ring_at_beginning_of_branch():
    """Test SELFIES that have a ring immediately at the start of a branch.
    """

    # CC1CCC(1CCl)F
    s = "[C][C][C][C][C][=Branch1][Branch1][Ring1][Ring2][C][Cl][F]"
    assert decode_eq(s, "CC1CCC1(CCl)F")

    # CC1CCS(Br)(1CCl)F
    s = "[C][C][C][C][S][Branch1][C][Br]" \
        "[=Branch1][Branch1][Ring1][Ring2][C][Cl][F]"
    assert decode_eq(s, "CC1CCS1(Br)(CCl)F")


def test_branch_and_ring_at_beginning_of_branch():
    """Test SELFIES that have a branch and ring immediately at the start
    of a branch.
    """

    # CC1CCCS((Br)1Cl)F
    s = "[C][C][C][C][C][S][#Branch1][#Branch1][Branch1][C][Br]" \
        "[Ring1][Branch1][Cl][F]"
    assert decode_eq(s, "CC1CCCS1(Br)(Cl)F")

    # CC1CCCS(1(Br)Cl)F
    s = "[C][C][C][C][C][S][#Branch1][#Branch1][Ring1][Branch1]" \
        "[Branch1][C][Br][Cl][F]"
    assert decode_eq(s, "CC1CCCS1(Br)(Cl)F")


def test_ring_immediately_following_branch():
    """Test SELFIES that have a ring immediately following after a branch.
    """

    # CCC1CCCC(OCO)1
    s = "[C][C][C][C][C][C][C][Branch1][Ring2][O][C][O][Ring1][Branch1]"
    assert decode_eq(s, "CCC1CCCC1OCO")

    # CCC1CCCC(OCO)(F)1
    s = "[C][C][C][C][C][C][C][Branch1][Ring2][O][C][O]" \
        "[Branch1][C][F][Ring1][Branch1]"
    assert decode_eq(s, "CCC1CCCC1(OCO)F")


def test_ring_after_branch():
    """Tests SELFIES that have a ring following a branch, but not
    immediately after a branch.
    """

    # CCCCCCC1(OCO)1
    s = "[C][C][C][C][C][C][C][Branch1][Ring2][O][C][O][C][Ring1][Branch1]"
    assert decode_eq(s, "CCCCCCC(OCO)=C")

    s = "[C][C][C][C][C][C][C][Branch1][Ring2][O][C][O]" \
        "[Branch1][C][F][C][C][Ring1][=Branch2]"
    assert decode_eq(s, "CCCCC1CC(OCO)(F)CC1")


def test_ring_on_top_of_existing_bond():
    """Tests SELFIES with rings between two atoms that are already bonded
    in the main scaffold.
    """

    # C1C1, C1C=1, C1C#1, ...
    assert decode_eq("[C][C][Ring1][C]", "C=C")
    assert decode_eq("[C][/C][Ring1][C]", "C=C")
    assert decode_eq("[C][C][=Ring1][C]", "C#C")
    assert decode_eq("[C][C][#Ring1][C]", "C#C")


def test_consecutive_rings():
    """Test SELFIES which have multiple consecutive rings.
    """

    s = "[C][C][C][C][Ring1][Ring2][Ring1][Ring2]"
    assert decode_eq(s, "C=1CCC=1")  # 1 + 1

    s = "[C][C][C][C][Ring1][Ring2][Ring1][Ring2][Ring1][Ring2]"
    assert decode_eq(s, "C#1CCC#1")  # 1 + 1 + 1

    s = "[C][C][C][C][=Ring1][Ring2][Ring1][Ring2]"
    assert decode_eq(s, "C#1CCC#1")  # 2 + 1

    s = "[C][C][C][C][Ring1][Ring2][=Ring1][Ring2]"
    assert decode_eq(s, "C#1CCC#1")  # 1 + 2

    # consecutive rings that exceed bond constraints
    s = "[C][C][C][C][#Ring1][Ring2][=Ring1][Ring2]"
    assert decode_eq(s, "C#1CCC#1")  # 3 + 2

    s = "[C][C][C][C][=Ring1][Ring2][#Ring1][Ring2]"
    assert decode_eq(s, "C#1CCC#1")  # 2 + 3

    s = "[C][C][C][C][=Ring1][Ring2][=Ring1][Ring2]"
    assert decode_eq(s, "C#1CCC#1")  # 2 + 2

    # consecutive rings with stereochemical single bond
    s = "[C][C][C][C][\\/Ring1][Ring2]"
    assert decode_eq(s, "C\\1CCC/1")

    s = "[C][C][C][C][\\/Ring1][Ring2][Ring1][Ring2]"
    assert decode_eq(s, "C=1CCC=1")


def test_unconstrained_symbols():
    """Tests SELFIES with symbols that are not semantically constrained.
    """

    f_branch = "[Branch1][C][F]"
    s = "[Xe-2]" + (f_branch * 8)
    assert decode_eq(s, "[Xe-2](F)(F)(F)(F)(F)(F)(F)CF")

    # change default semantic constraints
    constraints = sf.get_semantic_constraints()
    constraints["?"] = 2
    sf.set_semantic_constraints(constraints)

    assert decode_eq(s, "[Xe-2](F)CF")

    sf.set_semantic_constraints()


def test_isotope_symbols():
    """Tests that SELFIES symbols with isotope specifications are
     constrained properly.
    """

    s = "[13C][Branch1][C][Cl][Branch1][C][F][Branch1][C][Br][Branch1][C][I]"
    assert decode_eq(s, "[13C](Cl)(F)(Br)CI")

    assert decode_eq("[C][36Cl][C]", "C[36Cl]")


def test_chiral_symbols():
    """Tests that SELFIES symbols with chirality specifications are
    constrained properly.
    """

    s = "[C@@][Branch1][C][Cl][Branch1][C][F][Branch1][C][Br][Branch1][C][I]"
    assert decode_eq(s, "[C@@](Cl)(F)(Br)CI")

    s = "[C@H1][Branch1][C][Cl][Branch1][C][F][Branch1][C][Br]"
    assert decode_eq(s, "[C@H1](Cl)(F)CBr")


def test_explicit_hydrogen_symbols():
    """Tests that SELFIES symbols with explicit hydrogen specifications
     are constrained properly.
     """

    assert decode_eq("[CH1][Branch1][C][Cl][#C]", "[CH1](Cl)=C")
    assert decode_eq("[CH3][=C]", "[CH3]C")

    assert decode_eq("[CH4][C][C]", "[CH4]")
    assert decode_eq("[C][C][C][CH4]", "CCC")
    assert decode_eq("[C][Branch1][Ring2][C][=CH4][C][=C]", "C(C)=C")

    with pytest.raises(sf.DecoderError):
        sf.decoder("[C][C][CH5]")
    with pytest.raises(sf.DecoderError):
        sf.decoder("[C][C][C][OH9]")


def test_charged_symbols():
    """Tests that SELFIES symbols with charges are constrained properly.
    """

    constraints = sf.get_semantic_constraints()
    constraints["Sn+4"] = 1
    constraints["O-2"] = 2
    sf.set_semantic_constraints(constraints)

    # the following molecules don't make sense, but we use them to test
    # selfies. Hence, we can't verify them with RDKit
    assert decode_eq("[Sn+4][=C]", "[Sn+4]C")
    assert decode_eq("[O-2][#C]", "[O-2]=C")

    # mixing many symbol types
    assert decode_eq("[17O@@H1-2][#C]", "[17O@@H1-2]C")

    sf.set_semantic_constraints()


def test_standardized_alphabet():
    """Tests that equivalent SMILES atom symbols are translated into the
    same SELFIES atom symbol.
    """

    assert sf.encoder("[C][O][N][P][F]") == "[CH0][OH0][NH0][PH0][FH0]"
    assert sf.encoder("[Fe][Si]") == "[Fe][Si]"
    assert sf.encoder("[Fe++][Fe+2]") == "[Fe+2][Fe+2]"
    assert sf.encoder("[CH][CH1]") == "[CH1][CH1]"


def test_old_symbols():
    """Tests backward compatibility of SELFIES with old (<v2) symbols.
    """

    s = "[C@@Hexpl][Branch1_2][Branch1_1][Branch1_1][C][C][Cl][F]"
    assert sf.decoder(s, compatible=True) == "[C@@H1](C)(Cl)F"

    s = "[C][C][C][C][Expl=Ring1][Ring2][Expl#Ring1][Ring2]"
    assert sf.decoder(s, compatible=True) == "C#1CCC#1"

    long_s = "[C@@Hexpl][=C][C@@Hexpl][N+expl][=C][C+expl][N+expl][O+expl]" \
             "[Fe++expl][C@@Hexpl][C][N+expl][Branch1_2][Fe++expl][S+expl]" \
             "[=C][Expl=Ring1][Fe++expl][S+expl][Expl=Ring1][O+expl]" \
             "[C@@Hexpl][Expl=Ring1][C@@Hexpl][C@@Hexpl][N+expl][Expl=Ring1]" \
             "[Expl=Ring1][S+expl][=C]"
    try:
        sf.decoder(long_s, compatible=True)
    except Exception:
        assert False
