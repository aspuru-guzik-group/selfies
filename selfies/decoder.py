from collections import OrderedDict
from typing import Dict, Iterable, List, Optional, Tuple, Union

from selfies.grammar_rules import (get_bond_from_num,
                                   get_hypervalent_constraints,
                                   get_n_from_symbols, get_next_branch_state,
                                   get_next_state, get_num_from_bond,
                                   get_octet_rule_constraints,
                                   get_semantic_constraints,
                                   set_semantic_constraints)


def decoder(selfies: str,
            print_error: bool = False,
            constraints: Optional[str] = None) -> Optional[str]:
    """Translates a SELFIES into a SMILES.

    The SELFIES to SMILES translation operates based on the :mod:`selfies`
    grammar rules, which can be configured using
    :func:`selfies.set_semantic_constraints`. Given the appropriate settings,
    the decoded SMILES will always be syntactically and semantically correct.
    That is, the output SMILES will satisfy the specified bond constraints.
    Additionally, :func:`selfies.decoder` will attempt to preserve the
    atom and branch order of the input SELFIES.

    :param selfies: the SELFIES to be translated.
    :param print_error: if True, error messages will be printed to console.
        Defaults to False.
    :param constraints: if ``'octet_rule'`` or ``'hypervalent'``,
        the corresponding preset bond constraints will be used instead.
        If ``None``, :func:`selfies.decoder` will use the
        currently configured bond constraints. Defaults to ``None``.
    :return: the SMILES translation of ``selfies``. If an error occurs,
        and ``selfies`` cannot be translated, ``None`` is returned instead.

    :Example:

    >>> import selfies
    >>> selfies.decoder('[C][=C][F]')
    'C=CF'

    .. seealso:: The
        `"octet_rule" <https://en.wikipedia.org/wiki/Octet_rule>`_
        and
        `"hypervalent" <https://en.wikipedia.org/wiki/Hypervalent_molecule>`_
        preset bond constraints
        can be viewed with :func:`selfies.get_octet_rule_constraints` and
        :func:`selfies.get_hypervalent_constraints`, respectively. These
        presets are variants of the "default" bond constraints, which can
        be viewed with :func:`selfies.get_default_constraints`. Their
        differences can be summarized as follows:

            * def. : ``Cl``, ``Br``, ``I``: 1, ``N``: 3, ``P``: 5, ``P+1``: 6, ``P-1``: 4, ``S``: 6, ``S+1``: 7, ``S-1``: 5
            * oct. : ``Cl``, ``Br``, ``I``: 1, ``N``: 3, ``P``: 3, ``P+1``: 4, ``P-1``: 2, ``S``: 2, ``S+1``: 3, ``S-1``: 1
            * hyp. : ``Cl``, ``Br``, ``I``: 7, ``N``: 5, ``P``: 5, ``P+1``: 6, ``P-1``: 4, ``S``: 6, ``S+1``: 7, ``S-1``: 5
    """

    old_constraints = get_semantic_constraints()
    if constraints is None:
        pass
    elif constraints == 'octet_rule':
        set_semantic_constraints(get_octet_rule_constraints())
    elif constraints == 'hypervalent':
        set_semantic_constraints(get_hypervalent_constraints())
    else:
        raise ValueError("unrecognized constraint type")

    try:
        all_smiles = []  # process dot-separated fragments separately

        for s in selfies.split("."):
            smiles = _translate_selfies(s)

            if smiles != "":  # prevent malformed dots (e.g. [C]..[C], .[C][C])
                all_smiles.append(smiles)

        if constraints is not None:  # restore old constraints
            set_semantic_constraints(old_constraints)

        return '.'.join(all_smiles)

    except ValueError as err:
        if constraints is not None:  # restore old constraints
            set_semantic_constraints(old_constraints)

        if print_error:
            print("Decoding error '{}': {}.".format(selfies, err))
        return None


def _parse_selfies(selfies: str) -> Iterable[str]:
    """Parses a SELFIES into its symbols.

    A generator, which parses a SELFIES and yields its symbols
    one-by-one. When no symbols are left in the SELFIES, the empty
    string is infinitely yielded. As a precondition, the input SELFIES contains
    no dots, so all symbols are enclosed by square brackets, e.g. [X].

    :param selfies: the SElFIES string to be parsed.
    :return: an iterable of the symbols of the SELFIES.
    """

    left_idx = selfies.find('[')

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find(']', left_idx + 1)

        if (selfies[left_idx] != '[') or (right_idx == -1):
            raise ValueError("malformed SELIFES, "
                             "misplaced or missing brackets")

        next_symbol = selfies[left_idx: right_idx + 1]
        left_idx = right_idx + 1

        if next_symbol != '[nop]':  # skip [nop]
            yield next_symbol

    while True:  # no more symbols left
        yield ''


def _parse_selfies_symbols(selfies_symbols: List[str]) -> Iterable[str]:
    """Equivalent to ``_parse_selfies``, except the input SELFIES is presented
    as a list of SELFIES symbols, as opposed to a string.

    :param selfies_symbols: a SELFIES represented as a list of SELFIES symbols.
    :return: an iterable of the symbols of the SELFIES.
    """
    for symbol in selfies_symbols:

        if symbol != '[nop]':
            yield symbol

    while True:
        yield ''


def _translate_selfies(selfies: str) -> str:
    """A helper for ``selfies.decoder``, which translates a SELFIES into a
    SMILES (assuming the input SELFIES contains no dots).

    :param selfies: the SELFIES to be translated.
    :return: the SMILES translation of the SELFIES.
    """

    selfies_gen = _parse_selfies(selfies)

    # derived[i] is a list with three elements:
    #  (1) a string representing the i-th derived atom, and its connecting
    #      bond (e.g. =C, #N, N, C are all possible)
    #  (2) the number of available bonds the i-th atom has to make
    #  (3) the index of the previously derived atom that the i-th derived
    #      atom is bonded to
    # Example: if the 6-th derived atom was 'C', had 2 available bonds,
    # and was connected to the 5-th derived atom by a double bond, then
    # derived[6] = ['=C', 2, 5]
    derived = []

    # each item of <branches> is a key-value pair of indices that represents
    # the branches to be made. If a branch starts at the i-th derived atom
    # and ends at the j-th derived atom, then branches[i] = j. No two
    # branches should start at the same atom, e.g. C((C)Cl)C
    branches = {}

    # each element of <rings> is a tuple of size three that represents the
    # rings to be made, in the same order they appear in the SELFIES (left
    # to right). If the i-th ring is between the j-th and k-th derived atoms
    # (j <= k) and has bond symbol s ('=', '#', '\', etc.), then
    # rings[i] = (j, k, s).
    rings = []

    _translate_selfies_derive(selfies_gen, 0, derived, -1, branches, rings)
    _form_rings_bilocally(derived, rings)

    # create branches
    for lb, rb in branches.items():
        derived[lb][0] = '(' + derived[lb][0]
        derived[rb][0] += ')'

    smiles = ""
    for s, _, _ in derived:  # construct SMILES from <derived>
        smiles += s
    return smiles


# flake8: noqa: C901
# noinspection PyTypeChecker
def _translate_selfies_derive(selfies_gen: Iterable[str],
                              init_state: int,
                              derived: List[List[Union[str, int]]],
                              prev_idx: int,
                              branches: Dict[int, int],
                              rings: List[Tuple[int, int, str]]) -> None:
    """Recursive helper for _translate_selfies.

    Derives the SMILES symbols one-by-one from a SELFIES, and
    populates derived, branches, and rings. The main chain and side branches
    of the SELFIES are translated recursively. Rings are not actually
    translated, but saved to the rings list to be added later.

    :param selfies_gen: an iterable of the symbols of the SELFIES to be
        translated, created by ``_parse_selfies``.
    :param init_state: the initial derivation state.
    :param derived: see ``derived`` in ``_translate_selfies``.
    :param prev_idx: the index of the previously derived atom, or -1,
        if no atoms have been derived yet.
    :param branches: see ``branches`` in ``_translate_selfies``.
    :param rings: see ``rings`` in ``_translate_selfies``.
    :return: ``None``.
    """

    curr_symbol = next(selfies_gen)
    state = init_state

    while curr_symbol != '' and state >= 0:

        # Case 1: Branch symbol (e.g. [Branch1_2])
        if 'Branch' in curr_symbol:

            branch_init_state, new_state = \
                get_next_branch_state(curr_symbol, state)

            if state <= 1:  # state = 0, 1
                pass  # ignore no symbols

            else:
                L = int(curr_symbol[-4])  # corresponds to [BranchL_X]
                L_symbols = []
                for _ in range(L):
                    L_symbols.append(next(selfies_gen))

                N = get_n_from_symbols(*L_symbols)

                branch_symbols = []
                for _ in range(N + 1):
                    branch_symbols.append(next(selfies_gen))
                branch_gen = _parse_selfies_symbols(branch_symbols)

                branch_start = len(derived)
                _translate_selfies_derive(branch_gen, branch_init_state,
                                          derived, prev_idx, branches, rings)
                branch_end = len(derived) - 1

                # resolve C((C)Cl)C --> C(C)(Cl)C
                while branch_start in branches:
                    branch_start = branches[branch_start] + 1

                # finally, register the branch in branches
                if branch_start <= branch_end:
                    branches[branch_start] = branch_end

        # Case 2: Ring symbol (e.g. [Ring2])
        elif 'Ring' in curr_symbol:

            new_state = state

            if state == 0:
                pass  # ignore no symbols

            else:
                L = int(curr_symbol[-2])  # corresponds to [RingL]
                L_symbols = []
                for _ in range(L):
                    L_symbols.append(next(selfies_gen))

                N = get_n_from_symbols(*L_symbols)

                left_idx = max(0, prev_idx - (N + 1))
                right_idx = prev_idx

                bond_symbol = ''
                if curr_symbol[1:5] == 'Expl':
                    bond_symbol = curr_symbol[5]

                rings.append((left_idx, right_idx, bond_symbol))

        # Case 3: regular symbol (e.g. [N], [=C], [F])
        else:
            new_symbol, new_state = get_next_state(curr_symbol, state)

            if new_symbol != '':  # in case of [epsilon]
                derived.append([new_symbol, new_state, prev_idx])

                if prev_idx >= 0:
                    bond_num = get_num_from_bond(new_symbol[0])
                    derived[prev_idx][1] -= bond_num

                prev_idx = len(derived) - 1

        curr_symbol = next(selfies_gen)  # update symbol and state
        state = new_state


def _form_rings_bilocally(derived: List[List[Union[str, int]]],
                          rings: List[Tuple[int, int, str]]) -> None:
    """Forms all the rings specified by the rings list, in first-to-last order,
    by updating derived.

    :param derived: see ``derived`` in ``_translate_selfies``.
    :param rings: see ``rings`` in ``_translate_selfies``.
    :return: ``None``.
    """

    # due to the behaviour of allowing multiple rings between the same atom
    # pair, or rings between already bonded atoms, we first resolve all rings
    # so that only valid rings are left and placed into <ring_locs>.
    ring_locs = OrderedDict()

    for left_idx, right_idx, bond_symbol in rings:

        if left_idx == right_idx:  # ring to the same atom forbidden
            continue

        left_end = derived[left_idx]
        right_end = derived[right_idx]
        bond_num = get_num_from_bond(bond_symbol)

        if left_end[1] <= 0 or right_end[1] <= 0:
            continue  # no room for bond

        if bond_num > min(left_end[1], right_end[1]):
            bond_num = min(left_end[1], right_end[1])
            bond_symbol = get_bond_from_num(bond_num)

        # ring is formed between two atoms that are already bonded
        # e.g. CC1C1C --> CC=CC
        if left_idx == right_end[2]:

            right_symbol = right_end[0]

            if right_symbol[0] in {'-', '/', '\\', '=', '#'}:
                old_bond = right_symbol[0]
            else:
                old_bond = ''

            # update bond multiplicity and symbol
            new_bond_num = min(bond_num + get_num_from_bond(old_bond), 3)
            new_bond_symbol = get_bond_from_num(new_bond_num)

            right_end[0] = new_bond_symbol + right_end[0][len(old_bond):]

        # ring is formed between two atoms that are not bonded, e.g. C1CC1C
        else:
            loc = (left_idx, right_idx)

            if loc in ring_locs:
                # a ring is formed between two atoms that are have previously
                # been bonded by a ring, so ring bond multiplicity is updated

                new_bond_num = min(bond_num
                                   + get_num_from_bond(ring_locs[loc]), 3)
                new_bond_symbol = get_bond_from_num(new_bond_num)
                ring_locs[loc] = new_bond_symbol

            else:
                ring_locs[loc] = bond_symbol

        left_end[1] -= bond_num
        right_end[1] -= bond_num

    # finally, use <ring_locs> to add all the rings into <derived>

    ring_counter = 1
    for (left_idx, right_idx), bond_symbol in ring_locs.items():

        ring_id = str(ring_counter)
        if len(ring_id) == 2:
            ring_id = "%" + ring_id
        ring_counter += 1  # increment

        derived[left_idx][0] += bond_symbol + ring_id
        derived[right_idx][0] += bond_symbol + ring_id
