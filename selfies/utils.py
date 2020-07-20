from typing import Iterable, Set


def len_selfies(selfies: str) -> int:
    """Computes the character length of a SELFIES.

    The character length is the number of characters that make up the SELFIES,
    and not the length of the string itself (i.e. ``len(selfies)``).

    :param selfies: A SELFIES.
    :return: The character length of ``selfies``.

    :Example:

    >>> import selfies
    >>> selfies.len_selfies('[C][O][C]')
    3
    >>> selfies.len_selfies('[C][=C][F].[C]')
    5
    """

    return selfies.count('[') + selfies.count('.')


def split_selfies(selfies: str) -> Iterable[str]:
    """Splits a SELFIES into its characters.

    Returns an iterable that yields the characters of a SELFIES one-by-one
    in the order they appear in the string. SELFIES characters are always
    either indicated by an open and closed square bracket, or are the ``'.'``
    dot-bond character.

    :param selfies: The SELFIES to be read.
    :return: An iterable of the characters of ``selfies`` in the same order
        they appear in the string.

    :Example:

    >>> import selfies
    >>> list(selfies.split_selfies('[C][O][C]'))
    ['[C]', '[O]', '[C]']
    >>> list(selfies.split_selfies('[C][=C][F].[C]'))
    ['[C]', '[=C]', '[F]', '.', '[C]']
    """

    left_idx = selfies.find('[')

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find(']', left_idx + 1)
        next_char = selfies[left_idx: right_idx + 1]
        yield next_char

        left_idx = right_idx + 1
        if selfies[left_idx: left_idx + 1] == '.':
            yield '.'
            left_idx += 1


def get_alphabet_from_selfies(selfies_iter: Iterable[str]) -> Set[str]:
    """Constructs an alphabet from an iterable of SELFIES.

    From an iterable of SELFIES, constructs the minimum-sized set
    of SELFIES characters such that every SELFIES in the iterable can be
    constructed from characters from that set. Then, the set is returned.
    Note that the character ``.`` will not be added as a member of the
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
        for char in split_selfies(s):
            alphabet.add(char)

    alphabet.remove('.')

    return alphabet
