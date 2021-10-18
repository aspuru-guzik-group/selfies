import pytest

import selfies as sf


class Entry:

    def __init__(self, selfies, symbols, label, one_hot):
        self.selfies = selfies
        self.symbols = symbols
        self.label = label
        self.one_hot = one_hot


@pytest.fixture()
def dataset():
    stoi = {"[nop]": 0, "[O]": 1, ".": 2, "[C]": 3, "[F]": 4}
    itos = {i: c for c, i in stoi.items()}
    pad_to_len = 4

    entries = [
        Entry(selfies="",
              symbols=[],
              label=[0, 0, 0, 0],
              one_hot=[[1, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0]]),
        Entry(selfies="[C][C][C]",
              symbols=["[C]", "[C]", "[C]"],
              label=[3, 3, 3, 0],
              one_hot=[[0, 0, 0, 1, 0],
                       [0, 0, 0, 1, 0],
                       [0, 0, 0, 1, 0],
                       [1, 0, 0, 0, 0]]),
        Entry(selfies="[C].[C]",
              symbols=["[C]", ".", "[C]"],
              label=[3, 2, 3, 0],
              one_hot=[[0, 0, 0, 1, 0],
                       [0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0],
                       [1, 0, 0, 0, 0]]),
        Entry(selfies="[C][O][C][F]",
              symbols=["[C]", "[O]", "[C]", "[F]"],
              label=[3, 1, 3, 4],
              one_hot=[[0, 0, 0, 1, 0],
                       [0, 1, 0, 0, 0],
                       [0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 1]]),
        Entry(selfies="[C][O][C]",
              symbols=["[C]", "[O]", "[C]"],
              label=[3, 1, 3, 0],
              one_hot=[[0, 0, 0, 1, 0],
                       [0, 1, 0, 0, 0],
                       [0, 0, 0, 1, 0],
                       [1, 0, 0, 0, 0]])
    ]

    return entries, (stoi, itos, pad_to_len)


@pytest.fixture()
def dataset_flat_hots(dataset):
    flat_hots = []
    for entry in dataset[0]:
        hot = [elm for vec in entry.one_hot for elm in vec]
        flat_hots.append(hot)
    return flat_hots


def test_len_selfies(dataset):
    for entry in dataset[0]:
        assert sf.len_selfies(entry.selfies) == len(entry.symbols)


def test_split_selfies(dataset):
    for entry in dataset[0]:
        assert list(sf.split_selfies(entry.selfies)) == entry.symbols


def test_get_alphabet_from_selfies(dataset):
    entries, (vocab_stoi, _, _) = dataset

    selfies = [entry.selfies for entry in entries]
    alphabet = sf.get_alphabet_from_selfies(selfies)
    alphabet.add("[nop]")
    alphabet.add(".")

    assert alphabet == set(vocab_stoi.keys())


def test_selfies_to_encoding(dataset):
    entries, (vocab_stoi, vocab_itos, pad_to_len) = dataset

    for entry in entries:
        label, one_hot = sf.selfies_to_encoding(
            entry.selfies, vocab_stoi, pad_to_len, "both"
        )

        assert label == entry.label
        assert one_hot == entry.one_hot

        # recover original selfies
        selfies = sf.encoding_to_selfies(label, vocab_itos, "label")
        selfies = selfies.replace("[nop]", "")
        assert selfies == entry.selfies

        selfies = sf.encoding_to_selfies(one_hot, vocab_itos, "one_hot")
        selfies = selfies.replace("[nop]", "")
        assert selfies == entry.selfies


def test_selfies_to_flat_hot(dataset, dataset_flat_hots):
    entries, (vocab_stoi, vocab_itos, pad_to_len) = dataset

    batch = [entry.selfies for entry in entries]
    flat_hots = sf.batch_selfies_to_flat_hot(batch, vocab_stoi, pad_to_len)

    assert flat_hots == dataset_flat_hots

    # recover original selfies
    recovered = sf.batch_flat_hot_to_selfies(flat_hots, vocab_itos)
    assert batch == [s.replace("[nop]", "") for s in recovered]
