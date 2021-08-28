from selfies.exceptions import DecoderError
from selfies.grammar_rules import (
    get_index_from_selfies,
    next_atom_state,
    next_branch_state,
    next_ring_state,
    process_atom_symbol,
    process_branch_symbol,
    process_ring_symbol
)
from selfies.mol_graph import MolecularDFSTree
from selfies.utils.selfies_utils import split_selfies
from selfies.utils.smiles_utils import mol_to_smiles


def decoder(selfies: str) -> str:
    mol = MolecularDFSTree()

    rings = []
    for s in selfies.split("."):
        _derive_mol_from_symbols(
            symbol_iter=_tokenize_selfies(s),
            mol=mol,
            selfies=selfies,
            max_derive=float("inf"),
            init_state=0,
            root_atom=None,
            rings=rings
        )
    _form_rings_bilocally(mol, rings)
    return mol_to_smiles(mol)


def _tokenize_selfies(selfies):
    if isinstance(selfies, str):
        symbol_iter = split_selfies(selfies)
    elif isinstance(selfies, list):
        symbol_iter = selfies
    else:
        raise ValueError()  # should not happen

    try:
        for symbol in symbol_iter:
            if symbol == "[nop]":
                continue
            yield symbol
    except ValueError as err:
        raise DecoderError(str(err)) from None


def _derive_mol_from_symbols(
        symbol_iter, mol, selfies, max_derive,
        init_state, root_atom, rings
):
    n_derived = 0
    state = init_state
    prev_atom = root_atom

    while (state is not None) and (n_derived < max_derive):

        try:  # retrieve next symbol
            symbol = next(symbol_iter)
            n_derived += 1
        except StopIteration:
            break

        # Case 1: Branch symbol (e.g. [Branch1])
        if "ch" == symbol[-4:-2]:

            output = process_branch_symbol(symbol)
            if output is None:
                _raise_decoder_error(selfies, symbol)
            btype, n = output

            if state <= 1:
                next_state = state
            else:
                binit_state, next_state = next_branch_state(btype, state)

                Q = _read_index_from_selfies(symbol_iter, n_symbols=n)
                n_derived += n + _derive_mol_from_symbols(
                    symbol_iter, mol, selfies, (Q + 1),
                    init_state=binit_state, root_atom=prev_atom, rings=rings
                )

        # Case 2: Ring symbol (e.g. [Ring2])
        elif "ng" == symbol[-4:-2]:

            output = process_ring_symbol(symbol)
            if output is None:
                _raise_decoder_error(selfies, symbol)
            ring_type, n, stereo = output

            if state == 0:
                next_state = state
            else:
                ring_order, next_state = next_ring_state(ring_type, state)
                bond_info = (ring_order, stereo)

                Q = _read_index_from_selfies(symbol_iter, n_symbols=n)
                n_derived += n
                lidx = max(0, prev_atom.index - (Q + 1))
                rings.append((mol.get_atom(lidx), prev_atom, bond_info))

        # Case 3: [epsilon]
        elif "eps" in symbol:
            next_state = 0 if (state == 0) else None

        # Case 4: regular symbol (e.g. [N], [=C], [F])
        else:

            output = process_atom_symbol(symbol)
            if output is None:
                _raise_decoder_error(selfies, symbol)
            (bond_order, stereo), atom = output
            cap = atom.bonding_capacity

            bond_order, next_state = next_atom_state(bond_order, cap, state)
            if bond_order == 0:
                if state == 0:
                    mol.add_atom(atom, True)
            else:
                mol.add_atom(atom)
                src, dst = prev_atom.index, atom.index
                mol.add_bond(src=src, dst=dst, order=bond_order, stereo=stereo)
            prev_atom = atom

        if next_state is None:
            break
        state = next_state

    while n_derived < max_derive:  # consume remaining tokens
        try:
            next(symbol_iter)
            n_derived += 1
        except StopIteration:
            break

    return n_derived


def _raise_decoder_error(selfies, invalid_symbol):
    err_msg = "invalid symbol '{}'\n\tSELFIES: {}".format(
        invalid_symbol, selfies
    )
    raise DecoderError(err_msg)


def _read_index_from_selfies(symbol_iter, n_symbols):
    index_symbols = []
    for _ in range(n_symbols):
        try:
            index_symbols.append(next(symbol_iter))
        except StopIteration:
            index_symbols.append(None)
    return get_index_from_selfies(*index_symbols)


def _form_rings_bilocally(mol, rings):
    rings_made = [0] * len(mol)

    for latom, ratom, bond_info in rings:
        lidx, ridx = latom.index, ratom.index

        if lidx == ridx:  # ring to the same atom forbidden
            continue

        order, (lstereo, rstereo) = bond_info
        lfree = latom.bonding_capacity - mol.get_bond_count(lidx)
        rfree = ratom.bonding_capacity - mol.get_bond_count(ridx)

        if lfree <= 0 or rfree <= 0:
            continue  # no room for ring bond
        order = min(order, lfree, rfree)

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
