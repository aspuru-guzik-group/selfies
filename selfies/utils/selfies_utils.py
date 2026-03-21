from itertools import chain
from typing import Iterable, Iterator, Set


def len_selfies(selfies: str) -> int:
    """Returns the number of symbols in a given SELFIES string.

    :param selfies: a SELFIES string.
    :return: the symbol length of the SELFIES string.

    :Example:

    >>> import selfies as sf
    >>> sf.len_selfies("[C][=C][F].[C]")
    5
    """

    return selfies.count("[") + selfies.count(".")


def split_selfies(selfies: str) -> Iterator[str]:
    """Tokenizes a SELFIES string into its individual symbols.

    :param selfies: a SELFIES string.
    :return: the symbols of the SELFIES string one-by-one with order preserved.

    :Example:

    >>> import selfies as sf
    >>> list(sf.split_selfies("[C][=C][F].[C]"))
    ['[C]', '[=C]', '[F]', '.', '[C]']
    """
    splitted_on_dots = selfies.split(".")
    splitted_on_brackets = [s[1: -1].split("][") for s in splitted_on_dots]
    if "[" in splitted_on_brackets[-1][-1]:
        raise ValueError("malformed SELFIES string, hanging '[' bracket")
    symbols = ([f"[{x}]" for x in xs] for xs in splitted_on_brackets)
    interleaved = (
        x if i == 0 else ["."] + x for i, x in enumerate(symbols)
    )
    return chain.from_iterable(interleaved)


def get_alphabet_from_selfies(selfies_iter: Iterable[str]) -> Set[str]:
    """Constructs an alphabet from an iterable of SELFIES strings.

    The returned alphabet is the set of all symbols that appear in the
    SELFIES strings from the input iterable, minus the dot ``.`` symbol.

    :param selfies_iter: an iterable of SELFIES strings.
    :return: an alphabet of SELFIES symbols, built from the input iterable.

    :Example:

    >>> import selfies as sf
    >>> selfies_list = ["[C][F][O]", "[C].[O]", "[F][F]"]
    >>> alphabet = sf.get_alphabet_from_selfies(selfies_list)
    >>> sorted(list(alphabet))
    ['[C]', '[F]', '[O]']
    """

    alphabet = set()
    for s in selfies_iter:
        for symbol in split_selfies(s):
            alphabet.add(symbol)
    alphabet.discard(".")
    return alphabet
