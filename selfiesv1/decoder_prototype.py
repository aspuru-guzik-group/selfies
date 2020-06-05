from selfiesv1.state_dicts import \
    get_chars_index, get_next_branch_state, get_next_state


def decoder(selfies, N_restrict=True, print_error=True):
    """Converts a SELFIES string into its SMILES representation.

    Args:
        selfies: the SELFIES string to be decoded
        N_restrict: if True, nitrogen will be constrained to 3 bonds
        print_error: if True error messages will be printed to console

    Returns: the SMILES translation of <selfies>. If an error occurs, and
             <selfies> cannot be translated, -1 is returned instead.
    """

    if not isinstance(selfies, str):
        return -1

    try:
        all_selfies = []  # process dot-separated fragments separately
        for s in selfies.split("."):
            all_selfies.append(_translate_selfies(s, N_restrict))
        return '.'.join(all_selfies)

    except ValueError as err:
        if print_error:
            print(err)
            print("Could not decode SELFIES. Please contact authors.")
        return -1


def _parse_selfies(selfies):
    """A generator, which parses a SELFIES string and returns one-by-one
    its characters. We assume <selfies> has no dots in it, so a character
    is denoted by an open and closed square bracket, i.e. [X].

    Args:
        selfies: the SElFIES string to be parsed

    Returns: the characters of <selfies> one-by-one.
    """

    left_idx = selfies.find('[')

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find(']', left_idx + 1)
        next_char = selfies[left_idx: right_idx + 1]
        left_idx = right_idx + 1

        yield next_char

    while True:  # no more characters left
        yield ''


def _translate_selfies(selfies, N_restrict):
    """A helper for selfies.decoder, which converts a SELFIES string
    without dots into its SMILES representation.

    Args:
        selfies: the SELFIES string (without dots) to be decoded
        N_restrict: if True, nitrogen will be constrained to 3 bonds

    Returns: the SMILES translation of <selfies>.
    """

    # the i-th element of <derived> is list with three elements:
    #  (1) a string representing the i-th derived atom
    #  (2) the number of available bonds this atom has to make
    #  (3) a list of indexes j of the atoms the i-th derived atom is connected
    #      to, where j < i
    # For example, if the 6-th derived atom was 'C', had 2 available bonds,
    # and was connected to the 0 and 5-th derived atoms, the 6-th element would
    # be ['C', 2, [0, 5]]
    derived = []

    # each element of <branches> is a tuple of two indexes, indicating where
    # in <derived> a branch starts and ends (in this order).
    branches = []

    ring_counter = [1]  # a running counter of the ring numbers used

    _translate_selfies_derive(selfies, 0, N_restrict, derived, -1,
                              branches, ring_counter)

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


def _translate_selfies_derive(selfies, init_state, N_restrict,
                              derived, prev_idx, branches, counter):
    """Recursive helper for _translate_selfies. Derives the SMILES characters
    one-by-one from a SELFIES string, and fills <derived> and <branches>.

    Args:
        selfies: the SELFIES string (without dots) to be decoded
        init_state: the initial derivation state
        N_restrict: if True, nitrogen will be constrained to 3 bonds
        derived: see <derived> in _translate_selfies
        prev_idx: the index of the previously derived atom, or -1, if
                  no atoms have been derived yet
        branches: see <branches> in _translate_selfies
        counter: see <ring_counter> in _translate_selfies

    Returns: None.
    """

    selfies_gen = _parse_selfies(selfies)

    first_flag = True
    curr_char = next(selfies_gen)
    state = init_state

    while curr_char != '' and state >= 0:

        # Case 1: Branch character (e.g. [Branch1_2]
        if 'Branch' in curr_char:

            branch_init_state, new_state = \
                get_next_branch_state(curr_char, state)

            if state == 0 or state == 1:
                next(selfies_gen)  # ignore next characters

            elif state == 9991 or state == 9992 or state == 9993:
                pass  # ignore no characters

            else:
                L = int(curr_char[-4])  # corresponds to [BranchL_X]
                L_symbols = []
                for _ in range(L):
                    L_symbols.append(next(selfies_gen))

                N = get_chars_index(*L_symbols, default=1)

                branch_selfies = ""
                for _ in range(N):
                    branch_selfies += next(selfies_gen)

                derived[prev_idx][1] = new_state

                branch_start = len(derived)
                _translate_selfies_derive(branch_selfies, branch_init_state,
                                          N_restrict, derived, prev_idx,
                                          branches, counter)
                branch_end = len(derived) - 1

                new_state = derived[prev_idx][1]
                if new_state == 0:
                    new_state = -1

                if branch_start <= branch_end:
                    branches.append((branch_start, branch_end))

        # Case 2: Ring character (e.g. [Ring2])
        elif 'Ring' in curr_char:

            if state == 0:
                next(selfies_gen)  # ignore next character
                new_state = state

            elif state == 9991 or state == 9992 or state == 993:
                new_state = state  # ignore no characters

            else:
                L = int(curr_char[-2])  # corresponds to [RingL]
                L_symbols = []
                for _ in range(L):
                    L_symbols.append(next(selfies_gen))

                N = get_chars_index(*L_symbols, default=5)

                new_state = _form_ring(curr_char, state, derived, N, counter)

        # Case 3: regular character (e.g. [N], [=C], [F])
        else:
            new_char, new_state = get_next_state(curr_char, state, N_restrict)

            if new_char != '':  # exclude the case of [epsilon]
                derived.append([new_char, new_state, [prev_idx]])

                if not first_flag:
                    bond_num = _get_bond_num(new_char[0])
                    derived[prev_idx][1] -= bond_num
                else:
                    first_flag = False

                prev_idx = len(derived) - 1

        curr_char = next(selfies_gen)
        state = new_state


def _form_ring(ring_char, state, derived, N, ring_counter):
    left_idx = max(0, len(derived) - 1 - (N + 1))
    right_idx = len(derived) - 1

    if left_idx == right_idx or left_idx in derived[right_idx][2]:
        return state

    left_end = derived[left_idx]
    right_end = derived[right_idx]

    bond_char = ''
    bond_num = 1

    if ring_char[1:4] == 'Expl':
        bond_char = ring_char[5]
        bond_num = _get_bond_num(bond_char)

    if left_end[1] >= bond_num and right_end[1] >= bond_num:

        # form ring
        ring_id = str(ring_counter[0])
        if len(ring_id) == 2:
            ring_id = "%" + ring_id
        ring_id = bond_char + ring_id

        ring_counter[0] += 1  # increment

        left_end[0] += ring_id
        left_end[1] -= bond_num

        right_end[0] += ring_id
        right_end[1] -= bond_num

        derived[right_idx][2].append(left_idx)

        if state == bond_num:
            return -1
        return state - bond_num

    return state


def _get_bond_num(bond_char):
    """
    Gets the bond number of a SMILES representation of a bond.

    Args:
        bond_char: the bond character, e.g., '=', '-', '#'

    Returns: the number of bonds <bond_char> represents, or 1 if the
             bond character is unknown.
    """

    if bond_char == "=":
        return 2
    elif bond_char == "#":
        return 3
    else:
        return 1


if __name__ == '__main__':
    import selfies as sf

    mol = "CC1(CN1)C#N"

    print(sf.encoder(mol))
    print(decoder(sf.encoder(mol)))
    # print(sf.decoder(sf.encoder(mol)))
