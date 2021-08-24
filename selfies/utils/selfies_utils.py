from typing import Iterable, Iterator, Set


def len_selfies(selfies: str) -> int:
    """Returns the number of symbols in a given SELFIES.

    :param selfies:
    :return: the symbol length of ``selfies``.

    :Example:

    >>> import selfies as sf
    >>> sf.len_selfies("[C][=C][F].[C]")
    5
    """

    return selfies.count("[") + selfies.count(".")


def split_selfies(selfies: str) -> Iterator[str]:
    """Tokenizes a SELFIES into its individual symbols.

    :param selfies:
    :return: the symbols of ``selfies`` one-by-one and in the same order
        they appear in the string.

    :Example:

    >>> import selfies as sf
    >>> list(sf.split_selfies("[C][=C][F].[C]"))
    ['[C]', '[=C]', '[F]', '.', '[C]']
    """

    left_idx = selfies.find("[")

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find("]", left_idx + 1)
        if right_idx == -1:
            raise ValueError("malformed SELFIES, hanging '[' bracket")

        next_symbol = selfies[left_idx: right_idx + 1]
        yield next_symbol

        left_idx = right_idx + 1
        if selfies[left_idx: left_idx + 1] == ".":
            yield "."
            left_idx += 1


def get_alphabet_from_selfies(selfies_iter: Iterable[str]) -> Set[str]:
    """Constructs an alphabet from an iterable of SELFIES.

    The returned alphabet is the set of all symbols that appear in the
    SELFIES from the given iterable, minus the dot ``"."`` symbol.

    :param selfies_iter: an iterable of SELFIES.
    :return: an alphabet of SELFIES symbols, built from ``selfies_iter``.

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
