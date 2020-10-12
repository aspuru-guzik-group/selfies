import pytest

import selfies as sf


class Case:

    def __init__(self, selfies, length, symbols, label, one_hot):
        self.selfies = selfies
        self.length = length
        self.symbols = symbols
        self.label = label
        self.one_hot = one_hot


@pytest.fixture()
def test_cases():
    stoi = {"[nop]": 0, "[epsilon]": 1, ".": 2, "[C]": 3, "[F]": 4}
    itos = {i: c for c, i in stoi.items()}
    pad_to_len = 4

    cases = [
        Case(selfies="", length=0,
             symbols=[],
             label=[0, 0, 0, 0],
             one_hot=[[1, 0, 0, 0, 0],
                          [1, 0, 0, 0, 0],
                          [1, 0, 0, 0, 0],
                          [1, 0, 0, 0, 0]]),
        Case(selfies="[C][C][C]", length=3,
             symbols=["[C]", "[C]", "[C]"],
             label=[3, 3, 3, 0],
             one_hot=[[0, 0, 0, 1, 0],
                          [0, 0, 0, 1, 0],
                          [0, 0, 0, 1, 0],
                          [1, 0, 0, 0, 0]]),
        Case(selfies="[C].[C]", length=3,
             symbols=["[C]", ".", "[C]"],
             label=[3, 2, 3, 0],
             one_hot=[[0, 0, 0, 1, 0],
                          [0, 0, 1, 0, 0],
                          [0, 0, 0, 1, 0],
                          [1, 0, 0, 0, 0]]),
        Case(selfies="[C][epsilon][C][F]", length=4,
             symbols=["[C]", "[epsilon]", "[C]", "[F]"],
             label=[3, 1, 3, 4],
             one_hot=[[0, 0, 0, 1, 0],
                          [0, 1, 0, 0, 0],
                          [0, 0, 0, 1, 0],
                          [0, 0, 0, 0, 1]]),
        Case(selfies="[C][epsilon][C]", length=3,
             symbols=["[C]", "[epsilon]", "[C]"],
             label=[3, 1, 3, 0],
             one_hot=[[0, 0, 0, 1, 0],
                          [0, 1, 0, 0, 0],
                          [0, 0, 0, 1, 0],
                          [1, 0, 0, 0, 0]])
    ]

    return cases, (stoi, itos, pad_to_len)


@pytest.fixture()
def test_cases_flat_hots(test_cases):
    flat_hots = []
    for case in test_cases[0]:
        hot = [elm for vec in case.one_hot for elm in vec]
        flat_hots.append(hot)
    return flat_hots


def test_len_selfies(test_cases):
    for case in test_cases[0]:
        assert sf.len_selfies(case.selfies) == case.length


def test_split_selfies(test_cases):
    for case in test_cases[0]:
        assert list(sf.split_selfies(case.selfies)) == case.symbols


def test_get_alphabet_from_selfies(test_cases):
    case_list, (vocab_stoi, _, _) = test_cases

    selfies = [case.selfies for case in case_list]
    alphabet = sf.get_alphabet_from_selfies(selfies)
    alphabet.add("[nop]")
    alphabet.add(".")

    assert alphabet == set(vocab_stoi.keys())


def test_selfies_to_encoding(test_cases):
    case_list, (vocab_stoi, vocab_itos, pad_to_len) = test_cases

    for case in case_list:
        label, one_hot = sf.selfies_to_encoding(case.selfies, vocab_stoi,
                                                pad_to_len=pad_to_len,
                                                enc_type='both')
        assert label == case.label
        assert one_hot == case.one_hot

        # recover original selfies
        selfies = sf.encoding_to_selfies(label, vocab_itos,
                                         enc_type='label')
        selfies = selfies.replace("[nop]", "")
        assert selfies == case.selfies

        selfies = sf.encoding_to_selfies(one_hot, vocab_itos,
                                         enc_type='one_hot')
        selfies = selfies.replace("[nop]", "")
        assert selfies == case.selfies


def test_selfies_to_flat_hot(test_cases, test_cases_flat_hots):
    case_list, (vocab_stoi, vocab_itos, pad_to_len) = test_cases

    selfies_batch = [case.selfies for case in case_list]
    flat_hots = sf.batch_selfies_to_flat_hot(selfies_batch,
                                             vocab_stoi,
                                             pad_to_len)
    assert flat_hots == test_cases_flat_hots

    selfies_recover = sf.batch_flat_hot_to_selfies(flat_hots, vocab_itos)
    assert selfies_batch == [s.replace("[nop]", "") for s in selfies_recover]
