from typing import Dict, List, Optional, Tuple, Set

from selfies.state_library import build_state_dict, default_atom_dict


def get_alphabet() -> Set[str]:
    """Returns the alphabet that ``selfies`` is currently operating on.

    More specifically, the alphabet is the set of SELFIES characters that
    ``selfies`` recognizes and can apply semantic constraints to. ``selfies``
    initially operates upon a default alphabet, which can later be changed using
    ``selfies.set_alphabet``. After retrieving the alphabet, it is copied
    and returned as a set, i.e., mutating the returned set has no effect on
    the behaviour of ``selfies``.

    Although the characters ``'[epsilon]'`` and ``'.'`` are always
    recognized by ``selfies``, they will never be members of the returned set.

    :return: The alphabet that ``selfies`` is currently operating on.

    .. note:: In order to one-hot or integer encode SELFIES strings,
        ``selfies.get_alphabet`` may be a poor choice. This is because
        ``selfies`` often includes many characters in its alphabet that never
        actually appear in the input set. Instead,
        ``selfies.get_alphabet_from_selfies`` may be preferred.

    """

    global _state_library

    alphabet = set(_state_library[0].keys())
    alphabet.update([
        '[Branch1_1]', '[Branch1_2]', '[Branch1_3]', '[Ring1]',
        '[Branch2_1]', '[Branch2_2]', '[Branch2_3]', '[Ring2]',
        '[Branch3_1]', '[Branch3_2]', '[Branch3_3]', '[Ring3]',
    ])
    alphabet.remove('[?]')
    alphabet.remove('[epsilon]')
    alphabet.add('[nop]')

    return alphabet


def get_atom_dict() -> Dict[str, int]:
    """Returns the ``atom_dict`` that ``selfies`` is currently operating on.

    The ``atom_dict`` is the argument of the most recent call of
    ``selfies.set_alphabet``, or a default dictionary if the method has not
    been called yet. Once retrieved, it is copied and then returned. See
    ``selfies.set_alphabet`` for further explanation on ``atom_dict``.

    :return: The ``atom_dict`` that ``selfies`` is currently operating on.
    """

    global _atom_dict
    return dict(_atom_dict)


def set_alphabet(atom_dict: Optional[Dict[str, int]] = None) -> None:
    """Sets the alphabet the ``selfies`` is operating on based on **atom_dict**.

    The SELFIES alphabet and grammar is built dynamically from a dictionary
    **atom_dict** of atom(s) and/or ion(s) and their corresponding bond
    capacities. The key of the dictionary is a SMILES string representing
    either a single atom, or some atom(s) and/or ion(s) enclosed by square
    brackets. The corresponding value is the number of bonds that
    the key can make, between 1 and 8 inclusive. For example, one may have:

        * ``atom_dict['I'] = 1``
        * ``atom_dict['[C@@H]'] = 3``

    ``selfies.decoder`` will only generate SMILES that respect the bond
    constraints specified by the dictionary. In the example above, both
    ``'[C][=I]'`` and ``'[I][=C]'`` will be translated to ``'CI'`` and
    ``'IC'`` respectively, because ``I`` has been configured to make one bond
    maximally.

    If a SMILES key is not specified in **atom_dict**, it will by default be
    constrained to 8 bonds. To change the default setting for unrecognized
    characters, set ``atom_dict['?']`` to the desired integer (between 1 and 8
    inclusive). Note that ``selfies.decoder`` only recognizes the exact
    character string specified in **atom_dict**. For example, ``'[Fe+2]'`` will
    not be constrained if it is not in **atom_dict**, even if ``'[Fe++]'`` is
    a key in the dictionary.

    :param atom_dict: a dictionary of the atom(s) or ions that the new SELFIES
        alphabet will be built upon, with the value being the
        maximum bond capacity of the atom or ion. Defaults to ``None``. In that
        case, a default dictionary will be used for **atom_dict**.
    :return: ``None``.
    """

    global _state_library, _atom_dict
    _atom_dict = atom_dict
    _state_library = build_state_dict(atom_dict)


# Character State Dict Functions ===============================================

# <_state_library> is accessed through two keys, which are (1) the current
# derivation state and (2) the current SELFIES character to be derived, or
# '[?]' if the character is unknown. The corresponding value is a tuple of
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
        derived_char = _process_unknown_char(char)
        new_state = state_dict['[?]'][1] - get_num_from_bond(char[1])
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

_index_alphabet = ['[C]', '[Ring1]', '[Ring2]',
                   '[Branch1_1]', '[Branch1_2]', '[Branch1_3]',
                   '[Branch2_1]', '[Branch2_2]', '[Branch2_3]',
                    '[O]', '[N]', '[=N]', '[=C]', '[#C]', '[S]', '[P]']

# <_alphabet_code> takes as a key a SELFIES char, and its corresponding value
# is the index of the key.

_alphabet_code = {c: i for i, c in enumerate(_index_alphabet)}


def get_n_from_chars(*chars: List[str]) -> int:
    """Converts a list of SELFIES characters [c_1, ..., c_n] into a number N.
    This is done by converting each character c_n to an integer idx(c_n) via
    <_alphabet_code>, and then treating the list as a number in base
    len(_alphabet_code). If a character is unrecognized, it is given value 0 by
    default.

    Args:
        *chars: a list of SELFIES characters

    Returns: the corresponding N for <chars>
    """

    N = 0
    for i, c in enumerate(reversed(chars)):
        N_i = _alphabet_code.get(c, 0) * (len(_alphabet_code) ** i)
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
