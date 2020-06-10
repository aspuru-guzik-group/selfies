"""A file of various utilities, ranging from state dictionaries used to enforce
the SELFIES grammar in an efficient manner, helper methods, and the
selfies_alphabet() method.

Next steps include:
TODO: generate these _state_dicts and _state_library dynamically,
    using the above valence information.
TODO: For states 991-993, the new N state is 4, which is inconsistent with an
    unknown atom. Also, this can be expanded to pardon the restrictions on
    any atom in general.
TODO: add or remove error checking if needed
"""
from typing import List, Tuple


def selfies_alphabet() -> List[str]:
    """
    Returns: a list of the characters of the SELFIES alphabet
    """

    alphabet = ['[Branch1_1]', '[Branch1_2]', '[Branch1_3]', '[Ring1]',
                '[Branch2_1]', '[Branch2_2]', '[Branch2_3]', '[Ring2]',
                '[Branch3_1]', '[Branch3_2]', '[Branch3_3]', '[Ring3]',
                '[O]', '[=O]', '[S]', '[=S]',
                '[N]', '[=N]', '[#N]', '[NHexpl]', '[P]',
                '[C]', '[=C]', '[#C]',
                '[C@Hexpl]', '[C@@Hexpl]', '[C@expl]', '[C@@expl]',
                '[H]', '[F]']
    return alphabet


# Character State Dict Functions ===============================================

_state_dict_0 = {
    '[epsilon]': ('', 0),
    '[H]': ('[H]', 1),
    '[F]': ('F', 1),
    '[Cl]': ('Cl', 1),
    '[Br]': ('Br', 1),
    '[O]': ('O', 2),
    '[=O]': ('O', 2),
    '[N]': ('N', 3),
    '[=N]': ('N', 3),
    '[#N]': ('N', 3),
    '[NHexpl]': ('[NH]', 2),
    '[C]': ('C', 4),
    '[=C]': ('C', 4),
    '[#C]': ('C', 4),
    '[C@expl]': ('[C@]', 4),
    '[C@@expl]': ('[C@@]', 4),
    '[C@Hexpl]': ('[C@H]', 3),
    '[C@@Hexpl]': ('[C@@H]', 3),
    '[S]': ('S', 6),
    '[=S]': ('S', 6),
    '[???]': (None, 6)
}

_state_dict_1 = {
    '[epsilon]': ('', -1),
    '[H]': ('[H]', -1),
    '[F]': ('F', -1),
    '[Cl]': ('Cl', -1),
    '[Br]': ('Br', -1),
    '[O]': ('O', 1),
    '[=O]': ('O', -1),
    '[N]': ('N', 2),
    '[=N]': ('N', 2),
    '[#N]': ('N', 2),
    '[NHexpl]': ('[NH]', 1),
    '[C]': ('C', 3),
    '[=C]': ('C', 3),
    '[#C]': ('C', 3),
    '[C@expl]': ('[C@]', 3),
    '[C@@expl]': ('[C@@]', 3),
    '[C@Hexpl]': ('[C@H]', 2),
    '[C@@Hexpl]': ('[C@@H]', 2),
    '[S]': ('S', 5),
    '[=S]': ('S', 5),
    '[???]': (None, 6)
}

_state_dict_2 = {
    '[epsilon]': ('', -1),
    '[H]': ('[H]', -1),
    '[F]': ('F', -1),
    '[Cl]': ('Cl', -1),
    '[Br]': ('Br', -1),
    '[O]': ('O', 1),
    '[=O]': ('=O', -1),
    '[N]': ('N', 2),
    '[=N]': ('=N', 1),
    '[#N]': ('=N', 1),
    '[NHexpl]': ('[NH]', 1),
    '[C]': ('C', 3),
    '[=C]': ('=C', 2),
    '[#C]': ('=C', 2),
    '[C@expl]': ('[C@]', 3),
    '[C@@expl]': ('[C@@]', 3),
    '[C@Hexpl]': ('[C@H]', 2),
    '[C@@Hexpl]': ('[C@@H]', 2),
    '[S]': ('S', 5),
    '[=S]': ('=S', 4),
    '[???]': (None, 6)
}

_state_dict_3_to_6 = {
    '[epsilon]': ('', -1),
    '[H]': ('[H]', -1),
    '[F]': ('F', -1),
    '[Cl]': ('Cl', -1),
    '[Br]': ('Br', -1),
    '[O]': ('O', 1),
    '[=O]': ('=O', -1),
    '[N]': ('N', 2),
    '[=N]': ('=N', 1),
    '[#N]': ('#N', -1),
    '[NHexpl]': ('[NH]', 1),
    '[C]': ('C', 3),
    '[=C]': ('=C', 2),
    '[#C]': ('#C', 1),
    '[C@expl]': ('[C@]', 3),
    '[C@@expl]': ('[C@@]', 3),
    '[C@Hexpl]': ('[C@H]', 2),
    '[C@@Hexpl]': ('[C@@H]', 2),
    '[S]': ('S', 5),
    '[=S]': ('=S', 4),
    '[???]': (None, 6)
}

# <_state_library> is accessed through two keys, which are (1) the current
# derivation state and (2) the current SELFIES character to be derived, or
# '[???]' is the character is unknown. The corresponding value is a tuple of
# (1) the derived SMILES character, and (2) the next derivation state.

_state_library = {
    0: _state_dict_0,
    1: _state_dict_1,
    2: _state_dict_2,
    3: _state_dict_3_to_6,
    4: _state_dict_3_to_6,
    5: _state_dict_3_to_6,
    6: _state_dict_3_to_6,
    9991: _state_dict_1,
    9992: _state_dict_2,
    9993: _state_dict_3_to_6
}


def get_next_state(char: str, state: int, N_restrict: bool) -> Tuple[str, int]:
    """Given the current non-branch, non-ring character and current derivation
    state, retrieves the derived SMILES character and the next derivation state.

    Args:
        char: a SELFIES character that is not a Ring or Branch
        state: the current derivation state
        N_restrict: if True, nitrogen is restricted to 3 bonds

    Returns: a tuple of (1) the derived character, and (2) the
             next derivation state
    """

    state_dict = _state_library[state]

    if char in state_dict:
        derived_char, new_state = state_dict[char]

        # relax nitrogen constraints if N_restrict = False
        if not N_restrict and (char in ['[N]', '[=N]', ['#N']]):
            _, new_state = state_dict['[???]']

            if state >= 991:
                new_state = 4

        return derived_char, new_state

    else:  # unknown SELFIES character
        _, new_state = state_dict['[???]']
        derived_char = _process_unknown_char(char)
        return derived_char, new_state


# <_bracket_less_smiles> is a set of SELFIES symbols, whose
# SMILES counterparts cannot have brackets by convention.

_bracket_less_smiles = {'[B]', '[C]', '[N]', '[P]', '[O]', '[S]',
                        '[F]', '[Cl]', '[Br]', '[I]',
                        '[c]', '[n]', '[o]', '[s]', '[p]'}


def _process_unknown_char(char: str) -> str:
    """Attempts to convert an unknown SELFIES character <char> into a
    proper SMILES character. For example, explicit aromatic symbols
    are not part of the SELFIES alphabet, but _process_unknown_char
    will help convert [c], [n], [o] --> c, n, o.

    Args:
        char: an unknown SELFIES character

    Returns: the processed SMILES character
    """

    processed = ""

    if char[0: 2] in ('[=', '[#', '[\\', '[/', '[-'):
        processed += char[1]
        char = "[" + char[2:]

    if char in _bracket_less_smiles:
        char = char[1: -1]  # remove [ and ] brackets

    if 'expl' in char:
        char = char.replace('expl]', ']')

    processed += char

    return processed


# Branch State Dict Functions ==================================================

# <_branch_state_library> takes as a key the current derivation state.
# Its value is a tuple; for [BranchL_X], the (X - 1)th element of the tuple
# gives a tuple of (1) the initial branch derivation state and (2) the
# next derivation state (after the branch is derived). States 0-1, 9991-9993
# are not included because Branches at those states are simply skipped.

_branch_state_library = {
    0: ((-1, 0), (-1, 0), (-1, 0)),
    1: ((-1, 1), (-1, 1), (-1, 1)),
    2: ((9991, 1), (9991, 1), (9991, 1)),
    3: ((9991, 2), (9992, 1), (9992, 1)),
    4: ((9991, 3), (9992, 2), (9993, 1)),
    5: ((9991, 4), (9992, 3), (9993, 2)),
    6: ((9991, 5), (9992, 4), (9993, 3)),
    9991: ((-1, 9991), (-1, 9991), (-1, 9991)),
    9992: ((-1, 9992), (-1, 9992), (-1, 9992)),
    9993: ((-1, 9993), (-1, 9993), (-1, 9993))
}


def get_next_branch_state(branch_char: str, state: int) -> Tuple[int, int]:
    """Given the branch character and current derivation state, retrieves
    the initial branch derivation state (i.e. the derivation state that the
    new branch begins on), and the next derivation state (i.e. the derivation
    state after the branch is created).

    Args:
        branch_char: the branch character (e.g. [Branch1_2], [Branch3_1])
        state: the current derivation state

    Returns: a tuple of (1) the initial branch state and (2) the next state
    """

    branch_type = int(branch_char[-2])  # branches are of the form [BranchL_X]

    if not (1 <= branch_type <= 3):
        raise ValueError(f"Unknown branch character: {branch_char}")

    return _branch_state_library[state][branch_type - 1]


# SELFIES Character to N Functions =============================================

_index_alphabet = ['[epsilon]', '[Ring1]', '[Ring2]',
                   '[Branch1_1]', '[Branch1_2]', '[Branch1_3]',
                   '[Branch2_1]', '[Branch2_2]', '[Branch2_3]',
                   '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]',
                   '[C]', '[=C]', '[#C]', '[S]', '[=S]']

# <_alphabet_code> takes as a key a SELFIES char, and its corresponding value
# is the index of the key.

_alphabet_code = {c: i for i, c in enumerate(_index_alphabet)}


def get_n_from_chars(*chars: List[str], default: int = 1) -> int:
    """Converts a list of SELFIES characters [c_1, ..., c_n] into a number N.
    This is done by converting each character c_n to an integer idx(c_n) via
    <_alphabet_code>, and then treating the list as a number in base
    len(_alphabet_code).

    Args:
        *chars: a list of SELFIES characters
        default: the value to be returned if any character in <chars>
                 is not recognized

    Returns: the corresponding N for <chars>, or <default> if an element
             in <chars> does not have an index.
    """

    if any(c not in _alphabet_code for c in chars):
        return default

    N = 0
    for i, c in enumerate(reversed(chars)):
        N_i = _alphabet_code[c] * (len(_alphabet_code) ** i)
        N += N_i
    return N


def get_chars_from_n(n: int) -> List[str]:
    """Converts an integer n into a list of SELFIES characters that, if
    passed into <chars_to_index> in that order, would have produced n.

    Args:
        n: an integer

    Returns: a list of SELFIES characters representing n
             in base len(_alphabet_code)
    """

    if n == 0:
        return [_index_alphabet[0]]

    chars = []
    base = len(_index_alphabet)
    while n:
        chars.append(_index_alphabet[n % base])
        n //= base
    return chars[::-1]


# Helper Methods ===============================================================


def get_num_from_bond(bond_char: str) -> int:
    """Retrieves the bond multiplicity from a SMILES character representing
    a bond. If <bond_char> is not known, 1 is returned by default.

    Args:
        bond_char: a SMILES character representing a bond

    Returns: the bond multiplicity of <bond_char>, or 1 if <bond_char> is not
             recognized.
    """

    if bond_char == "=":
        return 2
    elif bond_char == "#":
        return 3
    else:
        return 1


def get_bond_from_num(n: int) -> str:
    """Returns the SMILES character representing a bond with multiplicity
    <n>. More specifically, '' = 1 and '=' = 2 and '#' = 3.

    Args:
        n: either 1, 2, 3

    Returns: the SMILES character representing a bond with multiplicity <n>
    """

    return ('', '=', '#')[n - 1]
