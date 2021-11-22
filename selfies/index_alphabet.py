import re
from typing import List, Dict, Union


_DEFAULT_INDEX_ALPHABET = {
    0: "[C]",
    1: "[Ring1]",
    2: "[Ring2]",
    3: "[Branch1]",
    4: "[=Branch1]",
    5: "[#Branch1]",
    6: "[Branch2]",
    7: "[=Branch2]",
    8: "[#Branch2]",
    9: "[O]",
    10: "[N]",
    11: "[=N]",
    12: "[=C]",
    13: "[#C]",
    14: "[S]",
    15: "[P]"
}

_PRESET_INDEX_ALPHABETS = {
    "default": dict(_DEFAULT_INDEX_ALPHABET),
}

_index_alphabet = _PRESET_INDEX_ALPHABETS["default"]
_index_alphabet_symbols = {
    _index_alphabet.get(c): c for i, c in enumerate(_index_alphabet)
}
_index_alphabet_reversed = tuple(
    [symbol for index_value, symbol in sorted(_index_alphabet.items())]
)


def get_preset_index_alphabet(name: str) -> Dict[int, str]:
    """Returns the preset index alphabet with the given name.
    :param name: the preset name: ``default`` or XXX or XXX.
    :return: the preset index alphabet with the specified name, represented as
        a dictionary which maps index values (the keys) to tokens (the values).
    """

    if name not in _PRESET_INDEX_ALPHABETS:
        raise ValueError("unrecognized preset name '{}'".format(name))
    return dict(_PRESET_INDEX_ALPHABETS[name])


def get_index_alphabet() -> Dict[int, str]:
    """Returns the semantic constraints that :mod:`selfies` is currently
    operating on.
    :return: the current semantic constraints, represented as a dictionary
        which maps index values (the keys) to tokens (the values).
    """

    global _index_alphabet
    return dict(_index_alphabet)


def update_index_alphabet(
        index_alphabet_updates: Union[str, Dict[int, str]] = "default"
) -> None:
    """Updates the index alphabet that :mod:`selfies` operates on.
    If the input is a string, the new index alphabet is taken to be
    the preset named ``index_alphabet_updates``
    (see :func:`selfies.get_preset_index_alphabet`).
    Otherwise, the input is a dictionary with updates for the index alphabet.
    This dictionary maps index values (the keys) to tokens (the values).
    :param index_alphabet_updates: the name of a preset, or a dictionary
        representing the new index alphabet.
    :return: ``None``.
    """

    global _index_alphabet
    global _index_alphabet_symbols
    global _index_alphabet_reversed

    SELFIES_ATOM_PATTERN = re.compile(
        r"^[\[]"  # opening square bracket [
        r"([=#/\\]?)"  # bond char
        r"(\d*)"  # isotope number (optional, e.g. 123, 26)
        r"([A-Z][a-z]?)"  # element symbol
        r"([@]{0,2})"  # chiral_tag (optional, only @ and @@ supported)
        r"((?:[H]\d)?)"  # H count (optional, e.g. H1, H3)
        r"((?:[+-][1-9]+)?)"  # charge (optional, e.g. +1)
        r"[]]$"  # closing square bracket ]
    )

    SELFIES_SPECIAL_TOKENS = set()
    for i in range(1, 4):
        SELFIES_SPECIAL_TOKENS.add("[Ring{}]".format(i))
        SELFIES_SPECIAL_TOKENS.add("[=Ring{}]".format(i))
        SELFIES_SPECIAL_TOKENS.add("[Branch{}]".format(i))
        SELFIES_SPECIAL_TOKENS.add("[=Branch{}]".format(i))
        SELFIES_SPECIAL_TOKENS.add("[#Branch{}]".format(i))

    if isinstance(index_alphabet_updates, str):
        _index_alphabet = get_preset_index_alphabet(index_alphabet_updates)
        _index_alphabet_reversed = {
            _index_alphabet.get(c): c for i, c in enumerate(_index_alphabet)
        }
        _index_alphabet = tuple(
            [symbol for index_value, symbol in sorted(_index_alphabet.items())]
        )

    elif isinstance(index_alphabet_updates, dict):

        _updated_index_alphabet = _index_alphabet.copy()
        _updated_index_alphabet.update(index_alphabet_updates)

        for key, value in _updated_index_alphabet.items():

            # error checking for index values
            if key not in [i for i in range(16)]:
                err_msg = "Invalid index value '{}' in index_alphabet_updates".format(key)
                raise ValueError(err_msg)

            # error checking for index symbols
            valid = False
            m = SELFIES_ATOM_PATTERN.match(value)
            if m is not None:
                valid = True
            if value in SELFIES_SPECIAL_TOKENS:
                valid = True
            if not valid:
                err_msg = "Invalid index symbol '{}' in index_alphabet_updates".format(value)
                raise ValueError(err_msg)

        # error checking for duplicate index symbols
        if not len(set(_updated_index_alphabet.values())) == 16:
            symbol_list = list(_updated_index_alphabet.values())
            err_msg = "Duplicate index symbol(s) '{}' in index_alphabet_updates".format(
                list(set([x for x in symbol_list if symbol_list.count(x) > 1]))
            )
            raise ValueError(err_msg)

        _index_alphabet = _updated_index_alphabet
        _index_alphabet_reversed = {
            _index_alphabet.get(c): c for i, c in enumerate(_index_alphabet)
        }
        _index_alphabet_symbols = tuple(
            [symbol for index_value, symbol in sorted(_index_alphabet.items())]
        )

    else:
        raise ValueError("index_alphabet must be a str or dict")


def get_index_from_selfies(*symbols: List[str]) -> int:
    index = 0
    base = len(_index_alphabet_reversed)
    for i, c in enumerate(reversed(symbols)):
        index += _index_alphabet_reversed.get(c, 0) * (base ** i)
    return index


def get_selfies_from_index(index: int) -> List[str]:
    if index < 0:
        raise IndexError()
    elif index == 0:
        return [_index_alphabet_symbols[0]]

    symbols = []
    base = len(_index_alphabet_symbols)
    while index:
        symbols.append(_index_alphabet_symbols[index % base])
        index //= base
    return symbols[::-1]
