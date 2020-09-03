from typing import Iterable, Set, Tuple, List, Union


def len_selfies(selfies: str) -> int:
    """Computes the symbol length of a SELFIES.

    The symbol length is the number of symbols that make up the SELFIES,
    and not the length of the string itself (i.e. ``len(selfies)``).

    :param selfies: A SELFIES.
    :return: The symbol length of ``selfies``.

    :Example:

    >>> import selfies
    >>> selfies.len_selfies('[C][O][C]')
    3
    >>> selfies.len_selfies('[C][=C][F].[C]')
    5
    """

    return selfies.count("[") + selfies.count(".")


def split_selfies(selfies: str) -> Iterable[str]:
    """Splits a SELFIES into its symbols.

    Returns an iterable that yields the symbols of a SELFIES one-by-one
    in the order they appear in the string. SELFIES symbols are always
    either indicated by an open and closed square bracket, or are the ``'.'``
    dot-bond symbol.

    :param selfies: The SELFIES to be read.
    :return: An iterable of the symbols of ``selfies`` in the same order
        they appear in the string.

    :Example:

    >>> import selfies
    >>> list(selfies.split_selfies('[C][O][C]'))
    ['[C]', '[O]', '[C]']
    >>> list(selfies.split_selfies('[C][=C][F].[C]'))
    ['[C]', '[=C]', '[F]', '.', '[C]']
    """

    left_idx = selfies.find("[")

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find("]", left_idx + 1)
        next_symbol = selfies[left_idx: right_idx + 1]
        yield next_symbol

        left_idx = right_idx + 1
        if selfies[left_idx: left_idx + 1] == ".":
            yield "."
            left_idx += 1


def get_alphabet_from_selfies(selfies_iter: Iterable[str]) -> Set[str]:
    """Constructs an alphabet from an iterable of SELFIES.

    From an iterable of SELFIES, constructs the minimum-sized set
    of SELFIES symbols such that every SELFIES in the iterable can be
    constructed from symbols from that set. Then, the set is returned.
    Note that the symbol ``'.'`` will not be added as a member of the
    returned set, even if it appears in the input.

    :param selfies_iter: An iterable of SELFIES.
    :return: The SElFIES alphabet built from the SELFIES in ``selfies_iter``.

    :Example:

    >>> import selfies
    >>> selfies_list = ['[C][F][O]', '[C].[O]', '[F][F]']
    >>> alphabet = selfies.get_alphabet_from_selfies(selfies_list)
    >>> sorted(list(alphabet))
    ['[C]', '[F]', '[O]']
    """

    alphabet = set()

    for s in selfies_iter:
        for symbol in split_selfies(s):
            alphabet.add(symbol)

    alphabet.discard(".")

    return alphabet


def selfies_to_hot(
    selfies: str, largest_selfie_len: int, alphabet: Union[List[str], Set[str]]
) -> Tuple[List[int], List[List[int]]]:
    """Go from a single selfies string to a one-hot encoding.

    :param selfies: SELFIES string to be converted.
    :param largest_selfie_len: Length of largest SELFIES used in
        construction of ``alphabet``.
    :param alphabet: Alphabet constructed from iterable of SELFIES.

    :return integer_encoded: padded, integer-encoded representation
        of ``selfies``.
    :return onehot_encoded: one-hot encoded representation with shape
        (largest_selfie_len, len(alphabet)).

    :Example:

    >>> import selfies as sf
    >>> sf.selfies_to_hot('[C]', 1, {'[C]'})
    ([0], [[1]])

    """

    # Should be a sorted list to preserve order of elements
    alphabet = sorted(list(alphabet))

    char_to_int = dict((c, i) for i, c in enumerate(alphabet))

    # pad with [nop]
    selfies += "[nop]" * (largest_selfie_len - len_selfies(selfies))

    # integer encode
    char_list = split_selfies(selfies)
    integer_encoded = [char_to_int[char] for char in char_list]

    # one-hot encode the integer encoded selfie
    onehot_encoded = list()
    for index in integer_encoded:
        letter = [0] * len(alphabet)
        letter[index] = 1
        onehot_encoded.append(letter)

    return (integer_encoded, onehot_encoded)


def multiple_selfies_to_hot(
    selfies_list: List[str],
    largest_selfie_len: int,
    alphabet: Union[List[str], Set[str]],
) -> List[List[int]]:
    """Convert a list of selfies strings to flattened one-hot encodings.

    :param selfies_list: SELFIES strings to be converted.
    :param largest_selfie_len: Length of largest SELFIES used in construction
        of ``alphabet``.
    :param alphabet: Alphabet constructed from iterable of SELFIES.

    :return hot_list: flattened one-hot encoded representations.

    :Example:

    >>> import selfies as sf
    >>> sl = ["[C]", "[C][C]"]
    >>> alpha = {'[C]', '[nop]'}
    >>> sf.multiple_selfies_to_hot(sl, 2, sorted(list(alpha)))
    [[1, 0, 0, 1], [1, 0, 1, 0]]

    """

    # Should be a sorted list to preserve order of elements
    alphabet = sorted(list(alphabet))

    hot_list = list()

    for selfies in selfies_list:
        _, onehot_encoded = selfies_to_hot(
            selfies, largest_selfie_len, alphabet
        )
        flattened = [elem for vec in onehot_encoded for elem in vec]
        hot_list.append(flattened)

    return hot_list


def hot_to_selfies(
    onehot_encoded: List[int],
    largest_selfie_len: int,
    alphabet: Union[List[str], Set[str]],
) -> str:
    """Convert a flattened one-hot encoding to a SELFIES string.

    :param onehot_encoded: one-hot encoded representation of a SELFIES.
    :param largest_selfie_len: Length of largest SELFIES used in construction
        of ``alphabet``.
    :param alphabet: Alphabet constructed from iterable of SELFIES.

    :return selfies: SELFIES string.

    :Example:

    >>> import selfies as sf
    >>> onehot = [1, 0, 0, 1, 0, 0, 1, 0, 0]
    >>> alpha = sorted(list({'[C]', '[epsilon]', '[nop]'}))
    >>> sf.hot_to_selfies(onehot, 3, alpha)
    '[C][C][C]'

    """

    # Should be a sorted list to preserve order of elements
    alphabet = sorted(list(alphabet))

    # Reshape to an N x M array where each column represents an alphabet
    # entry and each row is a position in the selfies
    twod_onehot = []
    N = largest_selfie_len
    M = len(alphabet)
    for i in range(N):
        twod_onehot.append(onehot_encoded[M * i: M * (i + 1)])

    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    integer_encoded = list()

    # Get integer encoding
    for row in twod_onehot:
        integer_encoded.append(row.index(1))

    # Integer encoding -> SELFIES
    char_list = [int_to_char[i] for i in integer_encoded]
    selfies = "".join(char_list)

    return selfies


def multiple_hot_to_selfies(
    onehot_encoded_list: List[List[int]],
    largest_selfie_len: int,
    alphabet: Union[List[str], Set[str]],
) -> List[str]:
    """Convert a list of one-hot encodings to SELFIES strings.

    :param onehot_encoded_list: one-hot encoded representations of SELFIES.
    :param largest_selfie_len: Length of largest SELFIES used in construction
        of ``alphabet``.
    :param alphabet: Alphabet constructed from iterable of SELFIES.

    :return selfies_list: SELFIES strings.

    :Example:

    >>> import selfies as sf
    >>> onehots = [[1, 0, 0, 1], [1, 0, 1, 0]]
    >>> alpha = {'[C]', '[nop]'}
    >>> sf.multiple_hot_to_selfies(onehots, 2, sorted(list(alpha)))
    ['[C][nop]', '[C][C]']

    """

    # Should be a sorted list to preserve order of elements
    alphabet = sorted(list(alphabet))

    selfies_list = []

    for onehot in onehot_encoded_list:
        selfies = hot_to_selfies(onehot, largest_selfie_len, alphabet)
        selfies_list.append(selfies)

    return selfies_list
