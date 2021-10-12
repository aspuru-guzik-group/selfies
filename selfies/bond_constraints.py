import functools
from itertools import product
from typing import Dict, Set, Union

from selfies.constants import ELEMENTS, INDEX_ALPHABET

_DEFAULT_CONSTRAINTS = {
    "H": 1, "F": 1, "Cl": 1, "Br": 1, "I": 1,
    "B": 3, "B+1": 2, "B-1": 4,
    "O": 2, "O+1": 3, "O-1": 1,
    "N": 3, "N+1": 4, "N-1": 2,
    "C": 4, "C+1": 5, "C-1": 3,
    "P": 5, "P+1": 6, "P-1": 4,
    "S": 6, "S+1": 7, "S-1": 5,
    "?": 8
}

_PRESET_CONSTRAINTS = {
    "default": dict(_DEFAULT_CONSTRAINTS),
    "octet_rule": dict(_DEFAULT_CONSTRAINTS),
    "hypervalent": dict(_DEFAULT_CONSTRAINTS)
}
_PRESET_CONSTRAINTS["octet_rule"].update(
    {"S": 2, "S+1": 3, "S-1": 1, "P": 3, "P+1": 4, "P-1": 2}
)
_PRESET_CONSTRAINTS["hypervalent"].update(
    {"Cl": 7, "Br": 7, "I": 7, "N": 5}
)

_current_constraints = _PRESET_CONSTRAINTS["default"]


def get_preset_constraints(name: str) -> Dict[str, int]:
    """Returns the preset semantic constraints with the given name.

    Besides the aforementioned default constraints, :mod:`selfies` offers
    other preset constraints for convenience; namely, constraints that
    enforce the `octet rule <https://en.wikipedia.org/wiki/Octet_rule>`_
    and constraints that accommodate `hypervalent molecules
    <https://en.wikipedia.org/wiki/Hypervalent_molecule>`_.

    The differences between these constraints can be summarized as follows:

    .. table::
        :align: center
        :widths: auto

        +-----------------+-----------+---+---+-----+-----+---+-----+-----+
        |                 | Cl, Br, I | N | P | P+1 | P-1 | S | S+1 | S-1 |
        +-----------------+-----------+---+---+-----+-----+---+-----+-----+
        | ``default``     |     1     | 3 | 5 |  6  |  4  | 6 |  7  |  5  |
        +-----------------+-----------+---+---+-----+-----+---+-----+-----+
        | ``octet_rule``  |     1     | 3 | 3 |  4  |  2  | 2 |  3  |  1  |
        +-----------------+-----------+---+---+-----+-----+---+-----+-----+
        | ``hypervalent`` |     7     | 5 | 5 |  6  |  4  | 6 |  7  |  5  |
        +-----------------+-----------+---+---+-----+-----+---+-----+-----+

    :param name: the preset name: ``default`` or ``octet_rule`` or
        ``hypervalent``.
    :return: the preset constraints with the specified name, represented
        as a dictionary which maps atoms (the keys) to their bonding capacities
        (the values).
    """

    if name not in _PRESET_CONSTRAINTS:
        raise ValueError("unrecognized preset name '{}'".format(name))
    return dict(_PRESET_CONSTRAINTS[name])


def get_semantic_constraints() -> Dict[str, int]:
    """Returns the semantic constraints that :mod:`selfies` is currently
    operating on.

    :return: the current semantic constraints, represented as a dictionary
        which maps atoms (the keys) to their bonding capacities (the values).
    """

    global _current_constraints
    return dict(_current_constraints)


def set_semantic_constraints(
        bond_constraints: Union[str, Dict[str, int]] = "default"
) -> None:
    """Updates the semantic constraints that :mod:`selfies` operates on.

    If the input is a string, the new constraints are taken to be
    the preset named ``bond_constraints``
    (see :func:`selfies.get_preset_constraints`).

    Otherwise, the input is a dictionary representing the new constraints.
    This dictionary maps atoms (the keys) to non-negative bonding
    capacities (the values); the atoms are specified by strings
    of the form ``E`` or ``E+C`` or ``E-C``,
    where ``E`` is an element symbol and ``C`` is a positive integer.
    For example, one may have:

       * ``bond_constraints["I-1"] = 0``
       * ``bond_constraints["C"] = 4``

    This dictionary must also contain the special ``?`` key, which indicates
    the bond capacities of all atoms that are not explicitly listed
    in the dictionary.

    :param bond_constraints: the name of a preset, or a dictionary
        representing the new semantic constraints.
    :return: ``None``.
    """

    global _current_constraints

    if isinstance(bond_constraints, str):
        _current_constraints = get_preset_constraints(bond_constraints)

    elif isinstance(bond_constraints, dict):

        # error checking
        if "?" not in bond_constraints:
            raise ValueError("bond_constraints missing '?' as a key")

        for key, value in bond_constraints.items():

            # error checking for keys
            j = max(key.find("+"), key.find("-"))
            if key == "?":
                valid = True
            elif j == -1:
                valid = (key in ELEMENTS)
            else:
                valid = (key[:j] in ELEMENTS) and key[j + 1:].isnumeric()
            if not valid:
                err_msg = "invalid key '{}' in bond_constraints".format(key)
                raise ValueError(err_msg)

            # error checking for values
            if not (isinstance(value, int) and value >= 0):
                err_msg = "invalid value at " \
                          "bond_constraints['{}'] = {}".format(key, value)
                raise ValueError(err_msg)

        _current_constraints = dict(bond_constraints)

    else:
        raise ValueError("bond_constraints must be a str or dict")

    # clear cache since we changed alphabet
    get_semantic_robust_alphabet.cache_clear()
    get_bonding_capacity.cache_clear()


@functools.lru_cache()
def get_semantic_robust_alphabet() -> Set[str]:
    """Returns a subset of all SELFIES symbols that are constrained
    by :mod:`selfies` under the current semantic constraints.

    :return: a subset of all SELFIES symbols that are semantically constrained.
    """

    alphabet_subset = set()
    bonds = {"": 1, "=": 2, "#": 3}

    # add atomic symbols
    for (a, c), (b, m) in product(_current_constraints.items(), bonds.items()):
        if (m > c) or (a == "?"):
            continue
        symbol = "[{}{}]".format(b, a)
        alphabet_subset.add(symbol)

    # add branch and ring symbols
    for i in range(1, 4):
        alphabet_subset.add("[Ring{}]".format(i))
        alphabet_subset.add("[=Ring{}]".format(i))
        alphabet_subset.add("[Branch{}]".format(i))
        alphabet_subset.add("[=Branch{}]".format(i))
        alphabet_subset.add("[#Branch{}]".format(i))

    alphabet_subset.update(INDEX_ALPHABET)

    return alphabet_subset


@functools.lru_cache()
def get_bonding_capacity(element: str, charge: int) -> int:
    """Returns the bonding capacity of a given atom, under the current
    semantic constraints.

    :param element: the element of the input atom.
    :param charge: the charge of the input atom.
    :return: the bonding capacity of the input atom.
    """

    key = element
    if charge != 0:
        key += "{:+}".format(charge)

    if key in _current_constraints:
        return _current_constraints[key]
    else:
        return _current_constraints["?"]
