from typing import Iterable, List


def len_selfies(selfies: str) -> int:
    """Retrieves the character length of a SELFIES. That is, the number
    of characters that make up the SELFIES; and not the length of the
    string itself (i.e. len(selfies))

    Args:
        selfies: a SELFIES

    Returns: the length of <selfies>
    """

    return selfies.count('[') + selfies.count('.')


def split_selfies(selfies: str) -> Iterable[str]:
    """Reads a SELFIES character-by-character.

    Args:
        selfies: a SELFIES

    Returns: an iterable of the characters of <selfies> in the same order they
             appear in the string
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


def get_alphabet_from_selfies(selfies_list: Iterable[str]) -> List[str]:
    """From a list of SELFIES, constructs the minimum-sized alphabet
    of SELFIES characters (e.g. '[C]', '.', '[Branch1_1]) such that
    every SELFIES in that list can be constructed from the alphabet.

    Args:
        selfies_list: a list of SELFIES

    Returns: a list of the characters of the SElFIES alphabet (in no
             particular oder)
    """

    alphabet = set()

    for s in selfies_list:
        for char in split_selfies(s):
            alphabet.add(char)

    return list(alphabet)
