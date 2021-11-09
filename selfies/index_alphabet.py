import functools
from itertools import product
from typing import Dict, Set, Union


_DEFAULT_INDEX_ALPHABET = {
    "[C]": 0, 
    "[Ring1]": 1, 
    "[Ring2]": 2, 
    "[Branch1]": 3, 
    "[=Branch1]": 4, 
    "[#Branch1]": 5, 
    "[Branch2]": 6, 
    "[=Branch2]": 7,
    "[#Branch2]": 8,
    "[O]": 9,
    "[N]": 10,
    "[=N]": 11,
    "[=C]": 12,
    "[#C]": 13,
    "[S]": 14,
    "[P]": 15
}

_PRESET_INDEX_ALPHABETS = {
    "default": dict(_DEFAULT_INDEX_ALPHABET),
}

_current_index_alphabet = _PRESET_INDEX_ALPHABETS["default"]

INDEX_ALPHABET = tuple(_current_index_alphabet)
INDEX_CODE = {c: i for i, c in enumerate(INDEX_ALPHABET)}

def get_preset_index_alphabets(name: str) -> Dict[str, int]:
    """the preset index alphabet with the given name.
    :param name: the preset name: ``default``.
    :return: the preset index alphabet with the specified name, represented as
        a dictionary which maps tokens (the keys) to their index values (the values).
    """

    if name not in _PRESET_INDEX_ALPHABETS:
        raise ValueError("unrecognized preset name '{}'".format(name))
    return dict(_PRESET_INDEX_ALPHABETS[name])

def get_current_index_alphabet() -> Dict[str, int]:
    """Returns the semantic constraints that :mod:`selfies` is currently
    operating on.
    :return: the current semantic constraints, represented as a dictionary
        which maps tokens (the keys) to their index values (the values).
    """

    global _current_index_alphabet
    return dict(_current_index_alphabet)


def set_index_alphabet(
        index_alphabet: Union[str, Dict[str, int]] = "default"
) -> None:
    """Updates the index alphabet that :mod:`selfies` operates on.
    If the input is a string, the new index alphabet is taken to be
    the preset named ``index_alphabet``
    (see :func:`selfies.get_preset_index_alphabets`).
    Otherwise, the input is a dictionary representing the new index alphabet.
    This dictionary maps tokens (the keys) to index values 
    (the values).
    :param index_alphabet: the name of a preset, or a dictionary
        representing the new index alphabet.
    :return: ``None``.
    """

    global _current_index_alphabet

    if isinstance(index_alphabet, str):
        _current_constraints = get_preset_index_alphabets(index_alphabet)

    elif isinstance(index_alphabet, dict):

        for key, value in index_alphabet.items():
            
            #NEED CODE TO DETERMINE IF SYMBOL IS VALID SELFIES SYMBOL
                
            if not (isinstance(value, int) and value >= 0 and value <= 16):
                err_msg = "invalid value at " \
                          "index_alphabet['{}'] = {}".format(key, value)
                raise ValueError(err_msg)
                

        _current_index_alphabet = dict(index_alphabet)

    else:
        raise ValueError("index_alphabet must be a str or dict")
