import pytest

import selfies as sf


@pytest.fixture()
def test_cases():
    return {
        "": (0, []),
        "[C][C][C]": (3, ["[C]", "[C]", "[C]"]),
        "[C].[C]": (3, ["[C]", ".", "[C]"]),
        "[C].[C][nop]": (4, ["[C]", ".", "[C]", "[nop]"]),
        "[C][epsilon][C]": (3, ["[C]", "[epsilon]", "[C]"]),
    }


@pytest.fixture()
def test_cases_alphabet():
    return {"[C]", "[nop]", "[epsilon]"}


@pytest.fixture()
def onehot_test_cases():
    return ["[C][C][C]", "[C][epsilon][C]", "[C][nop][nop]"]


@pytest.fixture()
def test_int_encodings():
    return [[0, 0, 0], [0, 1, 0], [0, 2, 2]]


@pytest.fixture()
def test_onehot_encodings():
    return [
        [1, 0, 0, 1, 0, 0, 1, 0, 0],
        [1, 0, 0, 0, 1, 0, 1, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0, 1],
    ]


def test_len_selfies(test_cases):
    for s, (length, _) in test_cases.items():
        assert sf.len_selfies(s) == length


def test_split_selfies(test_cases):
    for s, (_, symbols) in test_cases.items():
        assert list(sf.split_selfies(s)) == symbols


def test_get_alphabet_from_selfies(test_cases, test_cases_alphabet):
    alphabet = sf.get_alphabet_from_selfies(test_cases.keys())

    assert alphabet == test_cases_alphabet


def test_selfies_to_hot(
    onehot_test_cases, test_cases_alphabet, test_int_encodings
):
    alphabet = sorted(list(test_cases_alphabet))

    for idx, s in enumerate(onehot_test_cases):
        int_encoded, onehot_encoded = sf.selfies_to_hot(s, 3, alphabet)
        assert int_encoded == test_int_encodings[idx]
        assert len(onehot_encoded) == 3


def test_multiple_selfies_to_hot(
    onehot_test_cases, test_cases_alphabet, test_onehot_encodings
):
    alphabet = sorted(list(test_cases_alphabet))

    hot_list = sf.multiple_selfies_to_hot(onehot_test_cases, 3, alphabet)
    for v1, v2 in zip(test_onehot_encodings, hot_list):
        assert v1 == v2


def test_hot_to_selfies(
    onehot_test_cases, test_cases_alphabet, test_onehot_encodings
):
    alphabet = sorted(list(test_cases_alphabet))

    for idx, oh in enumerate(test_onehot_encodings):
        selfies = sf.hot_to_selfies(oh, 3, alphabet)
        assert selfies == onehot_test_cases[idx]


def test_multiple_hot_to_selfies(
    onehot_test_cases, test_cases_alphabet, test_onehot_encodings
):
    alphabet = sorted(list(test_cases_alphabet))

    selfies_list = sf.multiple_hot_to_selfies(
        test_onehot_encodings, 3, alphabet
    )
    assert onehot_test_cases == selfies_list
