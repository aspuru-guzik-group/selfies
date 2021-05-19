from rdkit import Chem

import selfies as sf


def reset_alphabet():
    sf.set_semantic_constraints({
        'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
        'O': 2, 'O+1': 3, 'O-1': 1,
        'N': 6, 'N+1': 4, 'N-1': 2,
        'C': 4, 'C+1': 5, 'C-1': 3,
        'S': 6, 'S+1': 7, 'S-1': 5,
        'P': 7, 'P+1': 8, 'P-1': 6,
        '?': 8,
    })


def test_branch_and_ring_at_state_X0():
    """Tests SELFIES with branches and rings at state X0 (i.e. at the
    very beginning of a SELFIES). These symbols should be skipped.
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[Branch3_1][C][S][C][O]"), "CSCO")
    assert is_eq(sf.decoder("[Ring3][C][S][C][O]"), "CSCO")
    assert is_eq(sf.decoder("[Branch1_1][Ring1][Ring3][C][S][C][O]"), "CSCO")


def test_branch_at_state_X1():
    """Test SELFIES with branches at state X1 (i.e. at an atom that
    can only make one bond. In this case, the branch symbol should be skipped.
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[C][C][O][Branch1_1][C][I]"), "CCOCI")
    assert is_eq(sf.decoder("[C][C][C][O][Branch3_3][C][I]"), "CCCOCI")


def test_branch_at_end_of_selfies():
    """Test SELFIES that have a branch symbol as its very last symbol.
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[C][C][C][C][Branch1_1]"), "CCCC")
    assert is_eq(sf.decoder("[C][C][C][C][Branch3_3]"), "CCCC")


def test_ring_at_end_of_selfies():
    """Test SELFIES that have a ring symbol as its very last symbol.
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[C][C][C][C][C][Ring1]"), "CCCC=C")
    assert is_eq(sf.decoder("[C][C][C][C][C][Ring3]"), "CCCC=C")


def test_branch_with_no_atoms():
    """Test SELFIES that have a branch, but the branch has no atoms in it.
    Such branches should not be made in the outputted SMILES.
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[C][Branch1_1][Ring2][Branch1_1]"
                            "[Branch1_1][Branch1_1][F]"),
                 "CF")
    assert is_eq(sf.decoder("[C][Branch1_1][Ring2][Ring1]"
                            "[Ring1][Branch1_1][F]"),
                 "CF")
    assert is_eq(sf.decoder("[C][Branch1_2][Ring2][Branch1_1]"
                            "[C][Cl][F]"),
                 "C(Cl)F")

    # special case: Branch3_3 takes Q_1, Q_2 = [O] and Q_3 = ''. However,
    # there are no more symbols in the branch.
    assert is_eq(sf.decoder("[C][C][C][C][Branch3_3][O][O]"), "CCCC")


def test_oversized_branch():
    """Test SELFIES that have a branch, with Q larger than the length
    of the SELFIES
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[C][Branch2_1][O][O][C][C][S][F][C]"), "C(CCSF)")
    assert is_eq(sf.decoder("[C][Branch2_3][O][O][#C][C][S][F]"), "C(#CCSF)")


def test_oversized_ring():
    """Test SELFIES that have a ring, with Q so large that the (Q + 1)-th
    previously derived atom does not exist.
    """

    reset_alphabet()
    assert is_eq(sf.decoder("[C][C][C][C][Ring1][O]"), "C1CCC1")
    assert is_eq(sf.decoder("[C][C][C][C][Ring2][O][C]"), "C1CCC1")

    # special case: Ring2 takes Q_1 = [O] and Q_2 = '', leading to
    # Q = 9 * 16 + 0 (i.e. an oversized ring)
    assert is_eq(sf.decoder("[C][C][C][C][Ring2][O]"), "C1CCC1")

    # special case: ring between 1st atom and 1st atom should not be formed
    assert is_eq(sf.decoder("[C][Ring1][O]"), "C")


def test_branch_at_beginning_of_branch():
    """Test SELFIES that have a branch immediately at the start of a branch.
    """

    reset_alphabet()

    # [C@]((Br)Cl)F
    assert is_eq(sf.decoder("[C@expl][Branch1_2][Branch1_1]"
                            "[Branch1_1][C][Br]"
                            "[Cl][F]"),
                 "[C@](Br)(Cl)F")

    # [C@](((Br)Cl)I)F
    assert is_eq(sf.decoder("[C@expl][Branch1_3][Branch2_1]"
                            "[Branch1_2][Branch1_1]"
                            "[Branch1_1][C][Br]"
                            "[Cl][I][F]"),
                 "[C@](Br)(Cl)(I)F")

    # [C@]((Br)(Cl)I)F
    assert is_eq(sf.decoder("[C@expl][Branch1_3][Branch2_1]"
                            "[Branch1_1][C][Br]"
                            "[Branch1_1][C][Cl]"
                            "[I][F]"),
                 "[C@](Br)(Cl)(I)F")


def test_ring_at_beginning_of_branch():
    """Test SELFIES that have a ring immediately at the start of a branch.
    """

    reset_alphabet()

    # CC1CCC(1CCl)F
    assert is_eq(sf.decoder("[C][C][C][C][C][Branch1_1][Branch1_1]"
                            "[Ring1][Ring2][C][Cl][F]"),
                 "CC1CCC1(CCl)F")

    # CC1CCS(Br)(1CCl)F
    assert is_eq(sf.decoder("[C][C][C][C][S][Branch1_1][C][Br]"
                            "[Branch1_1][Branch1_1][Ring1][Ring2][C][Cl][F]"),
                 "CC1CCS1(Br)(CCl)F")


def test_branch_and_ring_at_beginning_of_branch():
    """Test SELFIES that have a branch and ring immediately at the start
    of a branch.
    """

    reset_alphabet()

    # CC1CCCS((Br)1Cl)F
    assert is_eq(sf.decoder("[C][C][C][C][C][S][Branch1_2][Branch1_3]"
                            "[Branch1_1][C][Br]"
                            "[Ring1][Branch1_1][Cl][F]"),
                 "CC1CCCS1(Br)(Cl)F")

    # CC1CCCS(1(Br)Cl)F
    assert is_eq(sf.decoder("[C][C][C][C][C][S][Branch1_2][Branch1_3]"
                            "[Ring1][Branch1_1]"
                            "[Branch1_1][C][Br][Cl][F]"),
                 "CC1CCCS1(Br)(Cl)F")

    # CC1CCCS(((1Br)Cl)I)F
    assert is_eq(sf.decoder("[C][C][C][C][C][S][Branch1_3][Branch2_3]"
                            "[Branch1_2][Branch1_3]"
                            "[Branch1_1][Ring2][Ring1][Branch1_1][Br]"
                            "[Cl][I][F]"),
                 "CC1CCCS1(Br)(Cl)(I)F")


def test_ring_immediately_following_branch():
    """Test SELFIES that have a ring immediately following after a branch.
    """

    reset_alphabet()

    # CCC1CCCC(OCO)1
    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[Ring1][Branch1_1]"),
                 "CCC1CCCC1(OCO)")

    # CCC1CCCC(OCO)(F)1
    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[Branch1_1][C][F][Ring1][Branch1_1]"),
                 "CCC1CCCC1(OCO)(F)")


def test_ring_after_branch():
    """Tests SELFIES that have a ring following a branch, but not
    immediately after a branch.
    """

    reset_alphabet()

    # CCCCCCC1(OCO)1
    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[C][Ring1][Branch1_1]"),
                 "CCCCCCC(OCO)=C")

    assert is_eq(sf.decoder("[C][C][C][C][C][C][C][Branch1_1][Ring2][O][C][O]"
                            "[Branch1_1][C][F][C][C][Ring1][Branch2_2]"),
                 "CCCCC1CC(OCO)(F)CC1")


def test_ring_on_top_of_existing_bond():
    """Tests SELFIES with rings between two atoms that are already bonded
    in the main scaffold.
    """

    reset_alphabet()

    # C1C1, C1C=1, C1C#1, ...
    assert is_eq(sf.decoder("[C][C][Ring1][C]"), "C=C")
    assert is_eq(sf.decoder("[C][/C][Ring1][C]"), "C=C")
    assert is_eq(sf.decoder("[C][C][Expl=Ring1][C]"), "C#C")
    assert is_eq(sf.decoder("[C][C][Expl#Ring1][C]"), "C#C")


def test_consecutive_rings():
    """Test SELFIES which have multiple consecutive rings.
    """

    reset_alphabet()
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


def test_unconstrained_symbols():
    """Tests SELFIES with symbols that are not semantically constrained.
    """

    reset_alphabet()
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


def test_isotope_symbols():
    """Tests that SELFIES symbols with isotope specifications are
     constrained properly.
    """

    reset_alphabet()
    assert sf.decoder("[13Cexpl][Branch1_1][C][Cl][Branch1_1][C][F]"
                      "[Branch1_1][C][Br][Branch1_1][C][I]") \
           == "[13C](Cl)(F)(Br)CI"
    assert sf.decoder("[C][36Clexpl][C]") == "C[36Cl]"


def test_chiral_symbols():
    """Tests that SELFIES symbols with chirality specifications are
    constrained properly.
    """

    reset_alphabet()
    assert sf.decoder("[C@@expl][Branch1_1][C][Cl][Branch1_1][C][F]"
                      "[Branch1_1][C][Br][Branch1_1][C][I]") \
           == "[C@@](Cl)(F)(Br)CI"
    assert sf.decoder("[C@Hexpl][Branch1_1][C][Cl][Branch1_1][C][F]"
                      "[Branch1_1][C][Br]") \
           == "[C@H](Cl)(F)CBr"


def test_explicit_hydrogen_symbols():
    """Tests that SELFIES symbols with explicit hydrogen specifications
     are constrained properly.
     """

    reset_alphabet()
    assert sf.decoder("[CHexpl][Branch1_1][C][Cl][#C]") == "[CH](Cl)=C"
    assert sf.decoder("[CH3expl][=C]") == "[CH3]C"


def test_charged_symbols():
    """Tests that SELFIES symbols with charges are constrained properly.
    """

    reset_alphabet()
    constraints = sf.get_semantic_constraints()
    constraints['Sn+4'] = 1
    constraints['O-2'] = 2
    sf.set_semantic_constraints(constraints)

    # the following molecules don't make sense, but we use them to test
    # selfies. Hence, we can't verify them with RDKit
    assert sf.decoder("[Sn++++expl][=C]") == "[Sn++++]C"
    assert sf.decoder("[Sn+4expl][=C]") == "[Sn+4]C"
    assert sf.decoder("[O--expl][#C]") == "[O--]=C"
    assert sf.decoder("[O-2expl][#C]") == "[O-2]=C"

    # mixing many symbol types
    assert sf.decoder("[17O@@H-2expl][#C]") == "[17O@@H-2]C"

    sf.set_semantic_constraints()


# Helper Methods

def is_eq(smiles1, smiles2):
    return Chem.CanonSmiles(smiles1) == Chem.CanonSmiles(smiles2)


def test_decoder_infinite_loop_on_ring_order():
    """Tests that selfies containing interlocking loops, entered in non
     strictly increasing order, don't produce infinite looping.
    """
    smil = "O=C1NC6=C(C=CC=C6)C21OCC3(COC(C(NC5=C4C=CC=C5)=O)4OC3)CO2"
    equi = "O=C1NC6=C(C=CC=C6)C21OCC3(COC4(C(NC5=C4C=CC=C5)=O)OC3)CO2"
    result2 = sf.encoder(equi)
    result1 = sf.encoder(smil)
    assert result1 == result2
