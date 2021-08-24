import functools
from itertools import product
from typing import Dict, Set, Union

from selfies.constants import ELEMENTS, INDEX_ALPHABET

_DEFAULT_CONSTRAINTS = {
    "H": 1, "F": 1, "Cl": 1, "Br": 1, "I": 1,
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
    if name not in _PRESET_CONSTRAINTS:
        raise ValueError("unrecognized preset name '{}'".format(name))
    return dict(_PRESET_CONSTRAINTS[name])


def get_semantic_constraints() -> Dict[str, int]:
    global _current_constraints
    return dict(_current_constraints)


def set_semantic_constraints(
        bond_constraints: Union[str, Dict[str, int]] = "default"
) -> None:
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
            if not (isinstance(value, int) and value >= 1):
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
    key = element
    if charge != 0:
        key += "{:+}".format(charge)

    if key in _current_constraints:
        return _current_constraints[key]
    else:
        return _current_constraints["?"]
