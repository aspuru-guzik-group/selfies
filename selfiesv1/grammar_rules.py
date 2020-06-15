from typing import Dict, List, Optional, Tuple

from selfiesv1.state_library import build_state_dict, default_atom_dict


def get_alphabet() -> List[str]:
    """Retrieves the current SELFIES alphabet, which will be the default
    unless specified otherwise by <set_alphabet>.

    Returns: a list of the characters of the SELFIES alphabet (in no
             particular order)
    """

    global _state_library

    alphabet = list(_state_library[0].keys())
    alphabet.extend([
        '[Branch1_1]', '[Branch1_2]', '[Branch1_3]', '[Ring1]',
        '[Branch2_1]', '[Branch2_2]', '[Branch2_3]', '[Ring2]',
        '[Branch3_1]', '[Branch3_2]', '[Branch3_3]', '[Ring3]',
    ])
    alphabet.remove('[???]')

    return alphabet


def get_atom_dict() -> Dict[str, int]:
    """Retrieves the current <atom_dict> upon which the SELFIES alphabet is
    built upon. See <set_alphabet> for further explanation of
    the structure of <atom_dict>

    Returns: the current <atom_dict>
    """

    global _atom_dict
    return dict(_atom_dict)


def set_alphabet(atom_dict: Optional[Dict[str, int]] = None) -> None:
    """Sets the SELFIES alphabet to one based on the atom(s) or ions in
    <atom_dict>. <atom_dict> is a dictionary with the key being some atom(s)
    or ion represented in SMILES, and its corresponding value being the
    non-zero maximum bond capacity of the key. For example:
        atom_dict['C'] = 4
        atom_dict['Br'] = 1
        atom_dict['[C@@H]'] = 3
        atom_dict['[Cu++]'] = 4

    Args:
        atom_dict: a dictionary of the atoms or ions that the new SELFIES
                   alphabet will be built upon, with the value being the
                   maximum bond capacity of the atom or ion. If None,
                   then a default atom_dict will be used.

    Returns: None.
    """

    global _state_library, _atom_dict
    _atom_dict = atom_dict
    _state_library = build_state_dict(atom_dict)


# Character State Dict Functions ===============================================

# <_state_library> is accessed through two keys, which are (1) the current
# derivation state and (2) the current SELFIES character to be derived, or
# '[???]' is the character is unknown. The corresponding value is a tuple of
# (1) the derived SMILES character, and (2) the next derivation state.

_atom_dict = default_atom_dict
_state_library = build_state_dict()


def get_next_state(char: str, state: int) -> Tuple[str, int]:
    """Given the current non-branch, non-ring character and current derivation
    state, retrieves the derived SMILES character and the next derivation state.

    Args:
        char: a SELFIES character that is not a Ring or Branch
        state: the current derivation state

    Returns: a tuple of (1) the derived character, and (2) the
             next derivation state
    """

    state_dict = _state_library[state]

    if char in state_dict:
        return state_dict[char]

    else:  # unknown SELFIES character
        _, new_state = state_dict['[???]']
        derived_char = _process_unknown_char(char)
        return derived_char, new_state


# <_bracket_less_smiles> is a set of SELFIES symbols, whose
# SMILES counterparts cannot have brackets by convention.

_bracket_less_smiles = {'[B]', '[C]', '[N]', '[P]', '[O]', '[S]',
                        '[F]', '[Cl]', '[Br]', '[I]'}


def _process_unknown_char(char: str) -> str:
    """Attempts to convert an unknown SELFIES character <char> into a
    proper SMILES character.

    Args:
        char: an unknown SELFIES character

    Returns: the processed SMILES character
    """

    processed = ""

    if char[0: 2] in {'[=', '[#', '[\\', '[/', '[-'}:
        processed += char[1]
        char = "[" + char[2:]

    if char in _bracket_less_smiles:
        char = char[1: -1]  # remove [ and ] brackets

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
    7: ((9991, 6), (9992, 5), (9993, 4)),
    8: ((9991, 7), (9992, 6), (9993, 5)),
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

_index_alphabet = ['[C]', '[Ring1]', '[Ring2]', '[Ring3]',
                   '[Branch1_1]', '[Branch1_2]', '[Branch1_3]',
                   '[Branch2_1]', '[Branch2_2]', '[Branch2_3]',
                   '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]',
                   '[=C]', '[#C]', '[S]', '[=S]']

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
