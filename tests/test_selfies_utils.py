import pytest

import selfies as sf


@pytest.fixture()
def test_cases():
    return {
        "": (0, []),
        "[C][C][C]": (3, ["[C]", "[C]", "[C]"]),
        "[C].[C]": (3, ["[C]", ".", "[C]"]),
        "[C].[C][nop]": (4, ["[C]", ".", "[C]", "[nop]"]),
        "[C][epsilon][C]": (3, ["[C]", "[epsilon]", "[C]"])
    }


@pytest.fixture()
def test_cases_alphabet():
    return {'.', '[C]', '[nop]', '[epsilon]'}


def test_len_selfies(test_cases):
    for s, (length, _) in test_cases.items():
        assert sf.len_selfies(s) == length


def test_split_selfies(test_cases):
    for s, (_, symbols) in test_cases.items():
        assert list(sf.split_selfies(s)) == symbols


def test_get_alphabet_from_selfies(test_cases, test_cases_alphabet):
    alphabet = sf.get_alphabet_from_selfies(test_cases.keys())

    assert alphabet == test_cases_alphabet
