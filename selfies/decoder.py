from selfies.exceptions import DecoderError
from selfies.grammar_rules import (
    get_index_from_selfies,
    next_atom_state,
    next_branch_state,
    parse_atom_selfies,
    parse_branch_selfies,
    parse_ring_selfies
)
from selfies.mol_graph import MolecularDFSTree
from selfies.utils.selfies_utils import split_selfies
from selfies.utils.smiles_utils import mol_to_smiles


def decoder(selfies: str) -> str:
    smiles_fragments = []

    for s in selfies.split("."):
        mol = _parse_selfies(s)
        if len(mol) == 0:  # prevent malformed dots (e.g. [C]..[C], .[C][C])
            continue
        smiles_fragments.append(mol_to_smiles(mol))

    return ".".join(smiles_fragments)


def _tokenize_selfies(selfies):
    if isinstance(selfies, str):
        selfies_iter = split_selfies(selfies)
    elif isinstance(selfies, list):
        selfies_iter = selfies
    else:
        raise ValueError()  # should not happen

    try:
        for symbol in selfies_iter:
            if symbol == "[nop]":
                continue
            yield symbol
    except ValueError as err:
        raise DecoderError(str(err)) from None


def _parse_selfies(selfies):
    selfies_iter = _tokenize_selfies(selfies)
    mol = MolecularDFSTree()
    rings = []

    _derive_main_tree(
        mol, selfies, selfies_iter, float('inf'),
        init_state=-1, root_atom=None, rings=rings
    )

    _form_rings_bilocally(mol, rings)

    return mol


def _derive_main_tree(
        mol, selfies, selfies_iter, max_derive,
        init_state, root_atom, rings
):
    n_derived = 0
    state = init_state
    prev_atom = root_atom

    while (state != 0) and (n_derived < max_derive):

        try:  # retrieve next symbol
            symbol = next(selfies_iter)
            n_derived += 1
        except StopIteration:
            break

        # Case 1: Branch symbol (e.g. [Branch1])
        if 'ch' in symbol[-4:-2]:

            output = parse_branch_selfies(symbol)
            if output is None:
                _raise_decoder_error(selfies, symbol)
            branch_type, n = output
            br_init_state, new_state = next_branch_state(branch_type, state)

            if state > 1:
                Q = _read_index_from_selfies(selfies_iter, n_symbols=n)
                n_derived += n

                n_derived += _derive_main_tree(
                    mol, selfies, selfies_iter, (Q + 1),
                    init_state=br_init_state, root_atom=prev_atom, rings=rings
                )

        # Case 2: Ring symbol (e.g. [Ring2])
        elif 'ng' in symbol[-4:-2]:

            output = parse_ring_selfies(symbol)
            if output is None:
                _raise_decoder_error(selfies, symbol)
            bond_info, n = output
            new_state = state

            if state > 0:
                Q = _read_index_from_selfies(selfies_iter, n_symbols=n)
                n_derived += n

                lidx = max(0, prev_atom.index - (Q + 1))
                rings.append((mol.get_atom(lidx), prev_atom, bond_info))

        # Case 3: regular symbol (e.g. [N], [=C], [F])
        else:

            output = parse_atom_selfies(symbol)
            if output is None:
                _raise_decoder_error(selfies, symbol)
            (bond_order, stereo), atom = output
            cap = atom.bonding_capacity

            if (cap == 0) and (prev_atom is not None):
                _raise_decoder_error(selfies, symbol)
            bond_order, new_state = next_atom_state(bond_order, cap, state)

            mol.add_atom(atom)
            if prev_atom is not None:
                src, dst = prev_atom.index, atom.index
                mol.add_bond(src=src, dst=dst, order=bond_order, stereo=stereo)
            prev_atom = atom

        state = new_state

    while n_derived < max_derive:  # consume remaining tokens
        try:
            next(selfies_iter)
            n_derived += 1
        except StopIteration:
            break

    return n_derived


def _raise_decoder_error(selfies, invalid_symbol):
    err_msg = "invalid symbol '{}'\n\tSELFIES: {}".format(
        invalid_symbol, selfies
    )
    raise DecoderError(err_msg)


def _read_index_from_selfies(selfies_iter, n_symbols):
    index_symbols = []
    for _ in range(n_symbols):
        try:
            index_symbols.append(next(selfies_iter))
        except StopIteration:
            index_symbols.append(None)
    return get_index_from_selfies(*index_symbols)


def _form_rings_bilocally(mol, rings):
    rings_made = [0] * len(mol)

    for latom, ratom, bond_info in rings:
        lidx, ridx = latom.index, ratom.index

        if lidx == ridx:  # ring to the same atom forbidden
            continue

        (lorder, lstereo), (rorder, rstereo) = bond_info
        lfree = latom.bonding_capacity - mol.get_bond_count(lidx)
        rfree = ratom.bonding_capacity - mol.get_bond_count(ridx)

        if lfree <= 0 or rfree <= 0:
            continue  # no room for ring bond
        order = min(lorder, rorder, lfree, rfree)

        if mol.has_bond(a=lidx, b=ridx):
            bond = mol.get_dirbond(src=lidx, dst=ridx)
            new_order = min(order + bond.order, 3)
            mol.update_bond_order(a=lidx, b=ridx, new_order=new_order)

        else:
            mol.add_ring_bond(
                a=lidx, a_stereo=lstereo, a_pos=rings_made[lidx],
                b=ridx, b_stereo=rstereo, b_pos=rings_made[ridx],
                order=order
            )
            rings_made[lidx] += 1
            rings_made[ridx] += 1
