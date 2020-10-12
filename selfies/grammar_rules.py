from itertools import product
from typing import Dict, List, Optional, Set, Tuple

default_bond_constraints = {
    'H': 1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
    'O': 2, 'O+1': 3, 'O-1': 1,
    'N': 3, 'N+1': 4, 'N-1': 2,
    'C': 4, 'C+1': 5, 'C-1': 3,
    'S': 6, 'S+1': 7, 'S-1': 5,
    'P': 7, 'P+1': 8, 'P-1': 6,
    '?': 8,
}

_bond_constraints = default_bond_constraints


def get_semantic_robust_alphabet() -> Set[str]:
    """Returns a subset of all symbols that are semantically constrained
    by :mod:`selfies`.

    These semantic constraints can be configured with
    :func:`selfies.set_semantic_constraints`.

    :return: a subset of all symbols that are semantically constrained.
    """

    alphabet_subset = set()

    organic_subset = {'B', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'}
    bonds = {'': 1, '=': 2, '#': 3}

    # add atomic symbols
    for (a, c), (b, m) in product(_bond_constraints.items(), bonds.items()):

        if (m > c) or (a == '?'):
            continue

        if a in organic_subset:
            symbol = "[{}{}]".format(b, a)
        else:
            symbol = "[{}{}expl]".format(b, a)

        alphabet_subset.add(symbol)

    # add branch and ring symbols
    for i in range(1, 4):
        alphabet_subset.add("[Ring{}]".format(i))
        alphabet_subset.add("[Expl=Ring{}]".format(i))

        for j in range(1, 4):
            alphabet_subset.add("[Branch{}_{}]".format(i, j))

    return alphabet_subset


def get_semantic_constraints() -> Dict[str, int]:
    """Returns the semantic bond constraints that :mod:`selfies` is currently
    operating on.

    Returned is the argument of the most recent call of
    :func:`selfies.set_semantic_constraints`, or the default bond constraints
    if the function has not been called yet. Once retrieved, it is copied and
    then returned. See :func:`selfies.set_semantic_constraints` for further
    explanation.

    :return: the bond constraints :mod:`selfies` is currently operating on.
    """

    global _bond_constraints
    return dict(_bond_constraints)


def set_semantic_constraints(
        bond_constraints: Optional[Dict[str, int]] = None) -> None:
    """Configures the semantic constraints of :mod:`selfies`.

    The SELFIES grammar is enforced dynamically from a dictionary
    ``bond_constraints``. The keys of the dictionary are atoms and/or ions
    (e.g. ``I``, ``Fe+2``). To denote an ion, use the format ``E+C``
    or ``E-C``, where ``E`` is an element and ``C`` is a positive integer.
    The corresponding value is the maximum number of bonds that atom or
    ion can make, between 1 and 8 inclusive. For example, one may have:

        * ``bond_constraints['I'] = 1``
        * ``bond_constraints['C'] = 4``

    :func:`selfies.decoder` will only generate SMILES that respect the bond
    constraints specified by the dictionary. In the example above, both
    ``'[C][=I]'`` and ``'[I][=C]'`` will be translated to ``'CI'`` and
    ``'IC'`` respectively, because ``I`` has been configured to make one bond
    maximally.

    If an atom or ion is not specified in ``bond_constraints``, it will
    by default be constrained to 8 bonds. To change the default setting
    for unrecognized atoms or ions, set ``bond_constraints['?']`` to the
    desired integer (between 1 and 8 inclusive).

    :param bond_constraints: a dictionary representing the semantic
        constraints the updated SELFIES will operate upon. Defaults to
        ``None``; in this case, a default dictionary will be used.
    :return: ``None``.
    """

    global _bond_constraints

    if bond_constraints is None:
        _bond_constraints = default_bond_constraints

    else:

        # error checking
        if '?' not in bond_constraints:
            raise ValueError("bond_constraints missing '?' as a key.")

        for key, value in bond_constraints.items():
            if not (1 <= value <= 8):
                raise ValueError("bond_constraints['{}'] not between "
                                 "1 and 8 inclusive.".format(key))

        _bond_constraints = dict(bond_constraints)


# Symbol State Dict Functions ==============================================


def get_next_state(symbol: str, state: int) -> Tuple[str, int]:
    """Enforces the grammar rules for standard SELFIES symbols.

    Given the current non-branch, non-ring symbol and current derivation
    state, retrieves the derived SMILES symbol and the next derivation
    state.

    :param symbol: a SELFIES symbol that is not a Ring or Branch.
    :param state: the current derivation state.
    :return: a tuple of (1) the derived symbol, and
        (2) the next derivation state.
    """

    if symbol == '[epsilon]':
        return ('', 0) if state == 0 else ('', -1)

    # convert to smiles symbol
    bond = ''
    if symbol[1] in {'/', '\\', '=', '#'}:
        bond = symbol[1]
    bond_num = get_num_from_bond(bond)

    if symbol[-5:] == 'expl]':  # e.g. [C@@Hexpl]
        smiles_symbol = "[{}]".format(symbol[1 + len(bond):-5])
    else:
        smiles_symbol = symbol[1 + len(bond):-1]

    # get bond capacity
    element, h_count, charge = parse_atom_symbol(smiles_symbol)

    if charge == 0:
        atom_or_ion = element
    else:
        atom_or_ion = "{}{:+}".format(element, charge)

    max_bonds = _bond_constraints.get(atom_or_ion,
                                      _bond_constraints['?'])

    if (h_count > max_bonds) or (h_count == max_bonds and state > 0):
        raise ValueError("too many Hs in symbol '{}'; consider "
                         "adjusting bond constraints".format(symbol))
    max_bonds -= h_count  # hydrogens consume 1 bond

    # calculate next state
    if state == 0:
        bond = ''
        next_state = max_bonds
    else:
        if bond_num > min(state, max_bonds):
            bond_num = min(state, max_bonds)
            bond = get_bond_from_num(bond_num)

        next_state = max_bonds - bond_num
        if next_state == 0:
            next_state = -1

    return (bond + smiles_symbol), next_state


# Branch State Dict Functions =================================================


def get_next_branch_state(branch_symbol: str, state: int) -> Tuple[int, int]:
    """Enforces the grammar rules for SELFIES Branch symbols.

    Given the branch symbol and current derivation state, retrieves
    the initial branch derivation state (i.e. the derivation state that the
    new branch begins on), and the next derivation state (i.e. the derivation
    state after the branch is created).

    :param branch_symbol: the branch symbol (e.g. [Branch1_2], [Branch3_1])
    :param state: the current derivation state.
    :return: a tuple of (1) the initial branch state, and
        (2) the next derivation state.
    """

    branch_type = int(branch_symbol[-2])  # branches of the form [BranchL_X]

    if not (1 <= branch_type <= 3):
        raise ValueError("unknown branch symbol '{}'".format(branch_symbol))

    if 2 <= state <= 8:
        branch_init_state = min(state - 1, branch_type)
        next_state = state - branch_init_state
        return branch_init_state, next_state
    else:
        return -1, state


# SELFIES Symbol to N Functions ============================================

_index_alphabet = ['[C]', '[Ring1]', '[Ring2]',
                   '[Branch1_1]', '[Branch1_2]', '[Branch1_3]',
                   '[Branch2_1]', '[Branch2_2]', '[Branch2_3]',
                   '[O]', '[N]', '[=N]', '[=C]', '[#C]', '[S]', '[P]']

# _alphabet_code takes as a key a SELFIES symbol, and its corresponding value
# is the index of the key.

_alphabet_code = {c: i for i, c in enumerate(_index_alphabet)}


def get_n_from_symbols(*symbols: List[str]) -> int:
    """Computes N from a list of SELFIES symbols.

    Converts a list of SELFIES symbols [c_1, ..., c_n] into a number N.
    This is done by converting each symbol c_n to an integer idx(c_n) via
    ``_alphabet_code``, and then treating the list as a number in base
    len(_alphabet_code). If a symbol is unrecognized, it is given value 0 by
    default.

    :param symbols: a list of SELFIES symbols.
    :return: the corresponding N for ``symbols``.
    """

    N = 0
    for i, c in enumerate(reversed(symbols)):
        N_i = _alphabet_code.get(c, 0) * (len(_alphabet_code) ** i)
        N += N_i
    return N


def get_symbols_from_n(n: int) -> List[str]:
    """Converts an integer n into a list of SELFIES symbols that, if
    passed into ``get_n_from_symbols`` in that order, would have produced n.

    :param n: an integer from 0 to 4095 inclusive.
    :return: a list of SELFIES symbols representing n in base
        ``len(_alphabet_code)``.
    """

    if n == 0:
        return [_index_alphabet[0]]

    symbols = []
    base = len(_index_alphabet)
    while n:
        symbols.append(_index_alphabet[n % base])
        n //= base
    return symbols[::-1]


# Helper Functions ============================================================


def get_num_from_bond(bond_symbol: str) -> int:
    """Retrieves the bond multiplicity from a SMILES symbol representing
    a bond. If ``bond_symbol`` is not known, 1 is returned by default.

    :param bond_symbol: a SMILES symbol representing a bond.
    :return: the bond multiplicity of ``bond_symbol``, or 1 if
        ``bond_symbol`` is not recognized.
    """

    if bond_symbol == "=":
        return 2
    elif bond_symbol == "#":
        return 3
    else:
        return 1


def get_bond_from_num(n: int) -> str:
    """Returns the SMILES symbol representing a bond with multiplicity
    ``n``. More specifically, ``'' = 1`` and ``'=' = 2`` and ``'#' = 3``.

    :param n: either 1, 2, 3.
    :return: the SMILES symbol representing a bond with multiplicity ``n``.
    """

    return ('', '=', '#')[n - 1]


def find_element(atom_symbol: str) -> Tuple[int, int]:
    """Returns the indices of the element component of a SMILES atom symbol.

    That is, if atom_symbol[i:j] is the element substring of the SMILES atom,
    then (i, j) is returned. For example:
        *   _find_element('b') = (0, 1).
        *   _find_element('B') = (0, 1).
        *   _find_element('[13C]') = (3, 4).
        *   _find_element('[nH+]') = (1, 2).

    :param atom_symbol: a SMILES atom.
    :return: a tuple of the indices of the element substring of
        ``atom_symbol``.
    """

    if atom_symbol[0] != '[':
        return 0, len(atom_symbol)

    i = 1
    while atom_symbol[i].isdigit():  # skip isotope number
        i += 1

    if atom_symbol[i + 1].isalpha() and atom_symbol[i + 1] != 'H':
        return i, i + 2
    else:
        return i, i + 1


def parse_atom_symbol(atom_symbol: str) -> Tuple[str, int, int]:
    """Parses a SMILES atom symbol and returns its element component,
    number of associated hydrogens, and charge.

    See http://opensmiles.org/opensmiles.html for the formal grammar
    of SMILES atom symbols. Note that only @ and @@ are currently supported
    as chiral specifications.

    :param atom_symbol: a SMILES atom symbol.
    :return: a tuple of (1) the element of ``atom_symbol``, (2) the hydrogen
        count, and (3) the charge.
    """

    if atom_symbol[0] != '[':
        return atom_symbol, 0, 0

    atom_start, atom_end = find_element(atom_symbol)
    i = atom_end

    # skip chirality
    if atom_symbol[i] == '@':  # e.g. @
        i += 1
    if atom_symbol[i] == '@':  # e.g. @@
        i += 1

    h_count = 0  # hydrogen count
    if atom_symbol[i] == 'H':
        h_count = 1

        i += 1
        if atom_symbol[i].isdigit():  # e.g. [CH2]
            h_count = int(atom_symbol[i])
            i += 1

    charge = 0  # charge count
    if atom_symbol[i] in ('+', '-'):
        charge = 1 if atom_symbol[i] == '+' else -1

        i += 1
        if atom_symbol[i] in ('+', '-'):  # e.g. [Cu++]
            while atom_symbol[i] in ('+', '-'):
                charge += (1 if atom_symbol[i] == '+' else -1)
                i += 1

        elif atom_symbol[i].isdigit():  # e.g. [Cu+2]
            s = i
            while atom_symbol[i].isdigit():
                i += 1
            charge *= int(atom_symbol[s:i])

    return atom_symbol[atom_start: atom_end], h_count, charge
