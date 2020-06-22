from selfies.grammar_rules import get_bond_from_num, get_n_from_chars, \
    get_next_branch_state, get_next_state, get_num_from_bond

from typing import Optional, Iterable, List, Union, Tuple


def decoder(selfies: str, print_error: bool = False) -> Optional[str]:
    """Translates a SELFIES into a SMILES.

    The SELFIES to SMILES translation operates based on the ``selfies``
    alphabet and grammar, which can be changed using ``selfies.set_alphabet``.
    Given the appropriate settings, the decoded SMILES will always be
    syntactically and semantically correct. That is, the output SMILES will
    satisfy the specified bond constraints. Additionally, ``selfies.decoder``
    will attempt to preserve the atom and branch order of the input SELFIES.

    :param selfies: The SELFIES to be translated.
    :param print_error: If True, error messages will be printed to console.
        Defaults to False.
    :return: the SMILES translation of **selfies**. If an error occurs,
        and **selfies** cannot be translated, ``None`` is returned instead.

    :Example:

    >>> import selfies
    >>> selfies.decoder('[C][=C][F]')
    'C=CF'
    """

    try:
        all_smiles = []  # process dot-separated fragments separately

        for s in selfies.split("."):
            smiles = _translate_selfies(s)

            if smiles != "":  # prevent malformed dots (e.g. [C]..[C], .[C][C])
                all_smiles.append(smiles)

        return '.'.join(all_smiles)

    except ValueError as err:
        if print_error:
            print(err)
            print("Could not decode SELFIES. Please contact authors.")
        return None


def _parse_selfies(selfies: str) -> Iterable[str]:
    """Parses a SELFIES into its characters.

    A generator, which parses a SELFIES and yields its characters
    one-by-one. When no characters are left in the SELFIES, the empty
    string is infinitely yielded. As a precondition, the input SELFIES contains
    no dots, so all characters are enclosed by square brackets, e.g. [X].

    :param selfies: the SElFIES string to be parsed.
    :return: an iterable of the characters of the SELFIES.
    """

    left_idx = selfies.find('[')

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find(']', left_idx + 1)
        next_char = selfies[left_idx: right_idx + 1]
        left_idx = right_idx + 1

        if next_char != '[nop]':  # skip [nop]
            yield next_char

    while True:  # no more characters left
        yield ''


def _parse_selfies_chars(selfies_chars: List[str]) -> Iterable[str]:
    """Equivalent to ``_parse_selfies``, except the input SELFIES is presented
    as a list of SELFIES characters, as opposed to a string.

    :param selfies_chars: a list of SELFIES characters represented.
    :return: an iterable of the characters of the SELFIES.
    """
    for char in selfies_chars:

        if char != '[nop]':
            yield char

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

    # each element of <branches> is a tuple of size two that represents the
    # branches to be made, in the same order they appear in the SELFIES (left
    # to right). If the i-th branch starts at the j-th derived atom and ends
    # at the k-th derived atom, then branches[i] = (j, k).
    branches = []

    # each element of <rings> is a tuple of size three that represents the
    # rings to be made, in the same order they appear in the SELFIES (left
    # to right). If the i-th ring is between the j-th and k-th derived atoms
    # (j <= k) and has bond character s ('=', '#', '\', etc.), then
    # rings[i] = (j, k, s).
    rings = []

    _translate_selfies_derive(selfies_gen, 0, derived, -1, branches, rings)
    _form_rings_bilocally(derived, rings)

    lb_locs = {}  # key = index of left bracket, value = how many brackets
    rb_locs = {}  # key = index of right bracket, value = how many brackets
    for lb, rb in branches:
        lb_locs[lb] = lb_locs.get(lb, 0) + 1
        rb_locs[rb] = rb_locs.get(rb, 0) + 1

    smiles = ""
    for i, char in enumerate(derived):  # construct SMILES from <derived>

        if i in lb_locs:
            smiles += '(' * lb_locs[i]

        smiles += char[0]

        if i in rb_locs:
            smiles += ')' * rb_locs[i]

    return smiles


# flake8: noqa: C901
def _translate_selfies_derive(selfies_gen: Iterable[str],
                              init_state: int,
                              derived: List[List[Union[str, int]]],
                              prev_idx: int,
                              branches: List[Tuple[int, int]],
                              rings: List[Tuple[int, int, str]]) -> None:
    """Recursive helper for _translate_selfies.

    Derives the SMILES characters one-by-one from a SELFIES, and
    populates derived, branches, and rings. The main chain and side branches
    of the SELFIES are translated recursively. Rings are not actually
    translated, but saved to the rings list to be added later.

    :param selfies_gen: an iterable of the characters of the SELFIES to be
        translated, created by ``_parse_selfies``.
    :param init_state: the initial derivation state.
    :param derived: see ``derived`` in ``_translate_selfies``.
    :param prev_idx: the index of the previously derived atom, or -1,
        if no atoms have been derived yet.
    :param branches: see ``branches`` in ``_translate_selfies``.
    :param rings: see ``rings`` in ``_translate_selfies``.
    :return: ``None``.
    """

    curr_char = next(selfies_gen)
    state = init_state

    while curr_char != '' and state >= 0:

        # Case 1: Branch character (e.g. [Branch1_2])
        if 'Branch' in curr_char:

            branch_init_state, new_state = \
                get_next_branch_state(curr_char, state)

            if state <= 1 or state >= 9991:  # state = 0, 1, 9991, 9992, 9993
                pass  # ignore no characters

            else:
                L = int(curr_char[-4])  # corresponds to [BranchL_X]
                L_symbols = []
                for _ in range(L):
                    L_symbols.append(next(selfies_gen))

                N = get_n_from_chars(*L_symbols)

                branch_chars = []
                for _ in range(N + 1):
                    branch_chars.append(next(selfies_gen))
                branch_gen = _parse_selfies_chars(branch_chars)

                branch_start = len(derived)
                _translate_selfies_derive(branch_gen, branch_init_state,
                                          derived, prev_idx, branches, rings)
                branch_end = len(derived) - 1

                new_state = derived[prev_idx][1]
                if branch_start <= branch_end:
                    branches.append((branch_start, branch_end))

        # Case 2: Ring character (e.g. [Ring2])
        elif 'Ring' in curr_char:

            new_state = state

            if state == 0 or state >= 9991:  # state = 0, 9991, 9992, 9993
                pass  # ignore no characters

            else:
                L = int(curr_char[-2])  # corresponds to [RingL]
                L_symbols = []
                for _ in range(L):
                    L_symbols.append(next(selfies_gen))

                N = get_n_from_chars(*L_symbols)

                left_idx = max(0, len(derived) - 1 - (N + 1))
                right_idx = len(derived) - 1

                bond_char = ''
                if curr_char[1:5] == 'Expl':
                    bond_char = curr_char[5]

                rings.append((left_idx, right_idx, bond_char))

        # Case 3: regular character (e.g. [N], [=C], [F])
        else:
            new_char, new_state = get_next_state(curr_char, state)

            if new_char != '':  # in case of [epsilon]
                derived.append([new_char, new_state, prev_idx])

                if prev_idx >= 0:
                    bond_num = get_num_from_bond(new_char[0])
                    derived[prev_idx][1] -= bond_num

                prev_idx = len(derived) - 1

        curr_char = next(selfies_gen)  # update character and state
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
    ring_locs = {}

    for left_idx, right_idx, bond_char in rings:

        if left_idx == right_idx:  # ring to the same atom forbidden
            continue

        left_end = derived[left_idx]
        right_end = derived[right_idx]
        bond_num = get_num_from_bond(bond_char)

        if bond_num > left_end[1] or bond_num > right_end[1]:
            continue  # not enough available bonds to make the ring

        # ring is formed between two atoms that are already bonded
        # e.g. CC1C1C --> CC=CC
        if left_idx == right_end[2]:

            right_char = right_end[0]

            if right_char[0] in {'-', '/', '\\', '=', '#'}:
                old_bond = right_char[0]
            else:
                old_bond = ''

            # update bond multiplicity and character
            new_bond_num = min(bond_num + get_num_from_bond(old_bond), 3)
            new_bond_char = get_bond_from_num(new_bond_num)

            right_end[0] = new_bond_char + right_end[0][len(old_bond):]

        # ring is formed between two atoms that are not bonded, e.g. C1CC1C
        else:
            loc = (left_idx, right_idx)

            if loc in ring_locs:
                # a ring is formed between two atoms that are have previously
                # been bonded by a ring, so ring bond multiplicity is updated

                new_bond_num = min(bond_num
                                   + get_num_from_bond(ring_locs[loc]), 3)
                new_bond_char = get_bond_from_num(new_bond_num)
                ring_locs[loc] = new_bond_char

            else:
                ring_locs[loc] = bond_char

        left_end[1] -= bond_num
        right_end[1] -= bond_num

    # finally, use <ring_locs> to add all the rings into <derived>

    ring_counter = 1
    for (left_idx, right_idx), bond_char in ring_locs.items():

        ring_id = str(ring_counter)
        if len(ring_id) == 2:
            ring_id = "%" + ring_id
        ring_counter += 1  # increment

        derived[left_idx][0] += bond_char + ring_id
        derived[right_idx][0] += bond_char + ring_id
