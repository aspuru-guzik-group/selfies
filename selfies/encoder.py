from selfies.exceptions import EncoderError, SMILESParserError
from selfies.grammar_rules import get_selfies_from_index
from selfies.utils.smiles_utils import (
    atom_to_smiles,
    bond_to_smiles,
    smiles_to_mol
)

from selfies.mol_graph import AttributionMap


def encoder(smiles: str, strict: bool = True, attribute: bool = False) -> str:
    """Translates a SMILES string into its corresponding SELFIES string.

    This translation is deterministic and does not depend on the
    current semantic constraints. Additionally, it preserves the atom order
    of the input SMILES string; thus, one could generate randomized SELFIES
    strings by generating randomized SMILES strings, and then translating them.

    By nature of SELFIES, it is impossible to represent molecules that
    violate the current semantic constraints as SELFIES strings.
    Thus, we provide the ``strict`` flag to guard against such cases. If
    ``strict=True``, then this function will raise a
    :class:`selfies.EncoderError` if the input SMILES string represents
    a molecule that violates the semantic constraints. If
    ``strict=False``, then this function will not raise any error; however,
    calling :func:`selfies.decoder` on a SELFIES string generated this
    way will *not* be guaranteed to recover a SMILES string representing
    the original molecule.

    :param smiles: the SMILES string to be translated. It is recommended to
        use RDKit to check that the strings passed into this function
        are valid SMILES strings.
    :param strict: if ``True``, this function will check that the
        input SMILES string obeys the semantic constraints.
        Defaults to ``True``.
    :param attribute: if an attribution should be returned
    :return: a SELFIES string translated from the input SMILES string if
             attribute is ``False``, otherwise a tuple is returned of
             SELFIES string and attribution list.
    :raises EncoderError:  if the input SMILES string is invalid,
        cannot be kekulized, or violates the semantic constraints with
        ``strict=True``.

    :Example:

    >>> import selfies as sf
    >>> sf.encoder("C=CF")
    '[C][=C][F]'

    .. note:: This function does not currently support SMILES with:

        *   The wildcard symbol ``*``.
        *   The quadruple bond symbol ``$``.
        *   Chirality specifications other than ``@`` and ``@@``.
        *   Ring bonds across a dot symbol (e.g. ``c1cc([O-].[Na+])ccc1``) or
            ring bonds between atoms that are over 4000 atoms apart.

        Although SELFIES does not have aromatic symbols, this function
        *does* support aromatic SMILES strings by internally kekulizing them
        before translation.
    """

    try:
        mol = smiles_to_mol(smiles, attributable=attribute)
    except SMILESParserError as err:
        err_msg = "failed to parse input\n\tSMILES: {}".format(smiles)
        raise EncoderError(err_msg) from err

    if not mol.kekulize():
        err_msg = "kekulization failed\n\tSMILES: {}".format(smiles)
        raise EncoderError(err_msg)

    if strict:
        _check_bond_constraints(mol, smiles)

    # invert chirality of atoms where necessary,
    # such that they are restored when the SELFIES is decoded
    for atom in mol.get_atoms():
        if ((atom.chirality is not None)
                and mol.has_out_ring_bond(atom.index)
                and _should_invert_chirality(mol, atom)):
            atom.invert_chirality()

    fragments = []
    attribution_maps = []
    attribution_index = 0
    for root in mol.get_roots():
        derived = list(_fragment_to_selfies(
            mol, None, root, attribution_maps, attribution_index))
        attribution_index += len(derived)
        fragments.append("".join(derived))
    # trim attribution map of empty tokens
    attribution_maps = [a for a in attribution_maps if a.token]
    result = ".".join(fragments), attribution_maps
    return result if attribute else result[0]


def _check_bond_constraints(mol, smiles):
    errors = []

    for atom in mol.get_atoms():
        bond_cap = atom.bonding_capacity
        bond_count = mol.get_bond_count(atom.index)
        if bond_count > bond_cap:
            errors.append((atom_to_smiles(atom), bond_count, bond_cap))

    if errors:
        err_msg = "input violates the currently-set semantic constraints\n" \
                  "\tSMILES: {}\n" \
                  "\tErrors:\n".format(smiles)
        for e in errors:
            err_msg += "\t[{:} with {} bond(s) - " \
                       "a max. of {} bond(s) was specified]\n".format(*e)
        raise EncoderError(err_msg)


def _should_invert_chirality(mol, atom):
    out_bonds = mol.get_out_dirbonds(atom.index)

    # 1. rings whose right number are bonded to this atom (e.g. ...1...X1)
    # 2. rings whose left number are bonded to this atom (e.g. X1...1...)
    # 3. branches and other (e.g. X(...)...)
    partition = [[], [], []]
    for i, bond in enumerate(out_bonds):
        if not bond.ring_bond:
            partition[2].append(i)
        elif bond.src < bond.dst:
            partition[1].append(i)
        else:
            partition[0].append(i)
    partition[1].sort(key=lambda x: out_bonds[x].dst)

    # construct permutation
    perm = partition[0] + partition[1] + partition[2]
    count = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                count += 1
    return count % 2 != 0  # if odd permutation, should invert chirality


def _fragment_to_selfies(mol, bond_into_root, root,
                         attribution_maps, attribution_index=0):
    derived = []

    bond_into_curr, curr = bond_into_root, root
    while True:
        curr_atom = mol.get_atom(curr)
        token = _atom_to_selfies(bond_into_curr, curr_atom)
        derived.append(token)

        attribution_maps.append(AttributionMap(
            len(derived) - 1 + attribution_index,
            token, mol.get_attribution(curr_atom)))

        out_bonds = mol.get_out_dirbonds(curr)
        for i, bond in enumerate(out_bonds):

            if bond.ring_bond:
                if bond.src < bond.dst:
                    continue

                rev_bond = mol.get_dirbond(src=bond.dst, dst=bond.src)
                ring_len = bond.src - bond.dst
                Q_as_symbols = get_selfies_from_index(ring_len - 1)
                ring_symbol = "[{}Ring{}]".format(
                    _ring_bonds_to_selfies(rev_bond, bond),
                    len(Q_as_symbols)
                )

                derived.append(ring_symbol)
                attribution_maps.append(AttributionMap(
                    len(derived) - 1 + attribution_index,
                    ring_symbol, mol.get_attribution(bond)))
                for symbol in Q_as_symbols:
                    derived.append(symbol)
                    attribution_maps.append(AttributionMap(
                        len(derived) - 1 + attribution_index,
                        symbol, mol.get_attribution(bond)))

            elif i == len(out_bonds) - 1:
                bond_into_curr, curr = bond, bond.dst

            else:
                # start, end are so we can go back and
                # correct offset from branch symbol in
                # branch tokens
                start = len(attribution_maps)
                branch = _fragment_to_selfies(
                    mol, bond, bond.dst, attribution_maps, len(derived))
                Q_as_symbols = get_selfies_from_index(len(branch) - 1)
                branch_symbol = "[{}Branch{}]".format(
                    _bond_to_selfies(bond, show_stereo=False),
                    len(Q_as_symbols)
                )
                end = len(attribution_maps)

                derived.append(branch_symbol)
                for symbol in Q_as_symbols:
                    derived.append(symbol)
                    attribution_maps.append(AttributionMap(
                        len(derived) - 1 + attribution_index,
                        symbol, mol.get_attribution(bond)))

                # account for branch symbol because it is inserted after
                for j in range(start, end):
                    attribution_maps[j].index += len(Q_as_symbols) + 1
                attribution_maps.append(AttributionMap(
                    len(derived) - 1 + attribution_index,
                    branch_symbol, mol.get_attribution(bond)))

                derived.extend(branch)

        # end of chain
        if (not out_bonds) or out_bonds[-1].ring_bond:
            break
    return derived


def _bond_to_selfies(bond, show_stereo=True):
    if not show_stereo and (bond.order == 1):
        return ""
    return bond_to_smiles(bond)


def _ring_bonds_to_selfies(lbond, rbond):
    assert lbond.order == rbond.order

    if (lbond.order != 1) or all(b.stereo is None for b in (lbond, rbond)):
        return _bond_to_selfies(lbond, show_stereo=False)
    else:
        bond_char = "-" if (lbond.stereo is None) else lbond.stereo
        bond_char += "-" if (rbond.stereo is None) else rbond.stereo
        return bond_char


def _atom_to_selfies(bond, atom):
    assert not atom.is_aromatic
    bond_char = "" if (bond is None) else _bond_to_selfies(bond)
    return "[{}{}]".format(bond_char, atom_to_smiles(atom, brackets=False))
