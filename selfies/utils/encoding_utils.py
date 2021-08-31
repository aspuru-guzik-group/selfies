from typing import Dict, List, Tuple, Union

from selfies.utils.selfies_utils import len_selfies, split_selfies


def selfies_to_encoding(
        selfies: str,
        vocab_stoi: Dict[str, int],
        pad_to_len: int = -1,
        enc_type: str = 'both'
) -> Union[List[int], List[List[int]], Tuple[List[int], List[List[int]]]]:
    """Converts a SELFIES string into its label (integer)
    and/or one-hot encoding.

    A label encoded output will be a list of shape ``(L,)`` and a
    one-hot encoded output will be a 2D list of shape ``(L, len(vocab_stoi))``,
    where ``L`` is the symbol length of the SELFIES string. Optionally,
    the SELFIES string can be padded before it is encoded.

    :param selfies: the SELFIES string to be encoded.
    :param vocab_stoi: a dictionary that maps SELFIES symbols to indices,
        which must be non-negative and contiguous, starting from 0.
        If the SELFIES string is to be padded, then the special padding symbol
        ``[nop]`` must also be a key in this dictionary.
    :param pad_to_len: the length that the SELFIES string string is padded to.
        If this value is less than or equal to the symbol length of the
        SELFIES string, then no padding is added. Defaults to ``-1``.
    :param enc_type: the type of encoding of the output:
        ``label`` or ``one_hot`` or ``both``.
        If this value is ``both``, then a tuple of the label and one-hot
        encodings is returned. Defaults to ``both``.
    :return: the label encoded and/or one-hot encoded SELFIES string.

    :Example:

    >>> import selfies as sf
    >>> sf.selfies_to_encoding("[C][F]", {"[C]": 0, "[F]": 1})
    ([0, 1], [[1, 0], [0, 1]])
    """

    # some error checking
    if enc_type not in ("label", "one_hot", "both"):
        raise ValueError("enc_type must be in ('label', 'one_hot', 'both')")

    # pad with [nop]
    if pad_to_len > len_selfies(selfies):
        selfies += "[nop]" * (pad_to_len - len_selfies(selfies))

    # integer encode
    char_list = split_selfies(selfies)
    integer_encoded = [vocab_stoi[char] for char in char_list]

    if enc_type == "label":
        return integer_encoded

    # one-hot encode
    one_hot_encoded = list()
    for index in integer_encoded:
        letter = [0] * len(vocab_stoi)
        letter[index] = 1
        one_hot_encoded.append(letter)

    if enc_type == "one_hot":
        return one_hot_encoded
    return integer_encoded, one_hot_encoded


def encoding_to_selfies(
        encoding: Union[List[int], List[List[int]]],
        vocab_itos: Dict[int, str],
        enc_type: str,
) -> str:
    """Converts a label (integer) or one-hot encoding into a SELFIES string.

    If the input is label encoded, then a list of shape ``(L,)`` is
    expected; and if the input is one-hot encoded, then a 2D list of
    shape ``(L, len(vocab_itos))`` is expected.

    :param encoding: a label or one-hot encoding.
    :param vocab_itos: a dictionary that maps indices to SELFIES symbols.
        The indices of this dictionary must be non-negative and contiguous,
        starting from 0.
    :param enc_type: the type of encoding of the input:
        ``label`` or ``one_hot``.
    :return: the SELFIES string represented by the input encoding.

    :Example:

    >>> import selfies as sf
    >>> one_hot = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    >>> vocab_itos = {0: "[nop]", 1: "[C]", 2: "[F]"}
    >>> sf.encoding_to_selfies(one_hot, vocab_itos, enc_type="one_hot")
    '[C][F][nop]'
    """

    if enc_type not in ("label", "one_hot"):
        raise ValueError("enc_type must be in ('label', 'one_hot')")

    if enc_type == "one_hot":  # Get integer encoding
        integer_encoded = []
        for row in encoding:
            integer_encoded.append(row.index(1))
    else:
        integer_encoded = encoding

    # Integer encoding -> SELFIES
    char_list = [vocab_itos[i] for i in integer_encoded]
    selfies = "".join(char_list)

    return selfies


def batch_selfies_to_flat_hot(
        selfies_batch: List[str],
        vocab_stoi: Dict[str, int],
        pad_to_len: int = -1,
) -> List[List[int]]:
    """Converts a list of SELFIES strings into its list of flattened
    one-hot encodings.

    Each SELFIES string in the input list is one-hot encoded
    (and then flattened) using :func:`selfies.selfies_to_encoding`, with
    ``vocab_stoi`` and ``pad_to_len`` being passed in as arguments.

    :param selfies_batch: the list of SELFIES strings to be encoded.
    :param vocab_stoi: a dictionary that maps SELFIES symbols to indices.
    :param pad_to_len: the length that each SELFIES string in the input list
        is padded to. Defaults to ``-1``.
    :return: the flattened one-hot encodings of the input list.

    :Example:

    >>> import selfies as sf
    >>> batch = ["[C]", "[C][C]"]
    >>> vocab_stoi = {"[nop]": 0, "[C]": 1}
    >>> sf.batch_selfies_to_flat_hot(batch, vocab_stoi, 2)
    [[0, 1, 1, 0], [0, 1, 0, 1]]
    """

    hot_list = list()

    for selfies in selfies_batch:
        one_hot = selfies_to_encoding(selfies, vocab_stoi, pad_to_len,
                                      enc_type="one_hot")
        flattened = [elem for vec in one_hot for elem in vec]
        hot_list.append(flattened)

    return hot_list


def batch_flat_hot_to_selfies(
        one_hot_batch: List[List[int]],
        vocab_itos: Dict[int, str],
) -> List[str]:
    """Converts a list of flattened one-hot encodings into a list
    of SELFIES strings.

    Each encoding in the input list is unflattened and then decoded using
    :func:`selfies.encoding_to_selfies`, with ``vocab_itos`` being passed in
    as an argument.

    :param one_hot_batch: a list of flattened one-hot encodings. Each
        encoding must be a list of length divisible by ``len(vocab_itos)``.
    :param vocab_itos: a dictionary that maps indices to SELFIES symbols.
    :return: the list of SELFIES strings represented by the input encodings.

    :Example:

    >>> import selfies as sf
    >>> batch = [[0, 1, 1, 0], [0, 1, 0, 1]]
    >>> vocab_itos = {0: "[nop]", 1: "[C]"}
    >>> sf.batch_flat_hot_to_selfies(batch, vocab_itos)
    ['[C][nop]', '[C][C]']
    """

    selfies_list = []

    for flat_one_hot in one_hot_batch:

        # Reshape to an L x M array where each column represents an alphabet
        # entry and each row is a position in the selfies
        one_hot = []

        M = len(vocab_itos)
        if len(flat_one_hot) % M != 0:
            raise ValueError("size of vector in one_hot_batch not divisible "
                             "by the length of the vocabulary.")
        L = len(flat_one_hot) // M

        for i in range(L):
            one_hot.append(flat_one_hot[M * i: M * (i + 1)])

        selfies = encoding_to_selfies(one_hot, vocab_itos, enc_type="one_hot")
        selfies_list.append(selfies)

    return selfies_list
