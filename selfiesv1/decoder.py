from selfiesv1.utils import get_num_from_bond, get_n_from_chars, \
    get_next_branch_state, get_next_state, get_bond_from_num


def decoder(selfies, N_restrict=True, print_error=True):
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

    rings = []

    _translate_selfies_derive(selfies, 0, N_restrict, derived, -1,
                              branches, rings)
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


def _translate_selfies_derive(selfies, init_state, N_restrict,
                              derived, prev_idx, branches, rings):
    selfies_gen = _parse_selfies(selfies)

    curr_char = next(selfies_gen)
    state = init_state

    while curr_char != '' and state >= 0:

        # Case 1: Branch character (e.g. [Branch1_2])
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

                N = get_n_from_chars(*L_symbols, default=1)

                branch_selfies = ""
                for _ in range(N + 1):
                    branch_selfies += next(selfies_gen)

                branch_start = len(derived)
                _translate_selfies_derive(branch_selfies, branch_init_state,
                                          N_restrict, derived, prev_idx,
                                          branches, rings)
                branch_end = len(derived) - 1

                new_state = derived[prev_idx][1]
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

                N = get_n_from_chars(*L_symbols, default=5)

                left_idx = max(0, len(derived) - 1 - (N + 1))
                right_idx = len(derived) - 1

                bond_char = ''
                if curr_char[1:5] == 'Expl':
                    bond_char = curr_char[5]

                new_state = state
                rings.append((left_idx, right_idx, bond_char))

        # Case 3: regular character (e.g. [N], [=C], [F])
        else:
            new_char, new_state = get_next_state(curr_char, state, N_restrict)

            if new_char != '':  # exclude the case of [epsilon]
                derived.append([new_char, new_state, prev_idx])

                if prev_idx >= 0:
                    bond_num = get_num_from_bond(new_char[0])
                    derived[prev_idx][1] -= bond_num

                prev_idx = len(derived) - 1

        curr_char = next(selfies_gen)
        state = new_state


def _form_rings_bilocally(derived, rings):
    ring_locs = {}

    for left_idx, right_idx, bond_char in rings:

        if left_idx == right_idx:
            continue

        left_end = derived[left_idx]
        right_end = derived[right_idx]
        bond_num = get_num_from_bond(bond_char)

        if bond_num > left_end[1] or bond_num > right_end[1]:
            continue

        if left_idx == right_end[2]:

            right_char = right_end[0]

            if right_char[0] in {'-', '/', '\\', '=', '#'}:
                old_bond = right_char[0]
            else:
                old_bond = ''

            new_bond_num = min(bond_num + get_num_from_bond(old_bond), 3)
            new_bond_char = get_bond_from_num(new_bond_num)

            right_end[0] = new_bond_char + right_end[0][len(old_bond):]

        else:
            loc = (left_idx, right_idx)

            if loc in ring_locs:
                new_bond_num = min(bond_num
                                   + get_num_from_bond(ring_locs[loc]), 3)
                new_bond_char = get_bond_from_num(new_bond_num)
                ring_locs[loc] = new_bond_char

            else:
                ring_locs[loc] = bond_char

        left_end[1] -= bond_num
        right_end[1] -= bond_num

    ring_counter = 1
    for (left_idx, right_idx), bond_char in ring_locs.items():

        ring_id = str(ring_counter)
        if len(ring_id) == 2:
            ring_id = "%" + ring_id
        ring_counter += 1  # increment

        derived[left_idx][0] += bond_char + ring_id
        derived[right_idx][0] += bond_char + ring_id

