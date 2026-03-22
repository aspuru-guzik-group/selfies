from typing import Iterator
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


def _split_selfies_legacy(selfies: str) -> Iterator[str]:

    left_idx = selfies.find("[")

    while 0 <= left_idx < len(selfies):
        right_idx = selfies.find("]", left_idx + 1)
        if right_idx == -1:
            raise ValueError("malformed SELFIES string, hanging '[' bracket")

        next_symbol = selfies[left_idx: right_idx + 1]
        yield next_symbol

        left_idx = right_idx + 1
        if selfies[left_idx: left_idx + 1] == ".":
            yield "."
            left_idx += 1


def _vectorize_and_force_computation(func):
    def wrapper(xs):
        return [list(func(x)) for x in xs]
    return wrapper


class TestSplitSelfies:

    @pytest.fixture(scope="class")
    def datasets(self):
        smiles_list = """
            CC(=O)OCC[N+](C)(C)C
            CC(C[N+](C)(C)C)OC(=O)C
            O=C1CCCN1CC#CC[N+](C)(C)C
            NC(=O)OCC[N+](C)(C)C
            CC(C[N+](C)(C)C)OC(=O)N
            Cc1ccc(o1)C[N+](C)(C)C
            COC(=O)C1=CCCN(C1)C
            O=C1CCCN1CC#CCN1CCCC1
            CON=CC1=CCCN(C1)C
            CN1CC(=CCC1)C(=O)OCC#C
            CCC1C(=O)OCC1Cc1cncn1C
            CC(=O)OC1C[NH+]2CCC1CC2
        """.split()
        selfies_list = [sf.encoder(smiles) for smiles in smiles_list]
        selfies_with_dots = [f"{s}.{s}" for s in selfies_list]
        long_selfies = [s * 10 for s in selfies_list]
        return {
            "regular": selfies_list,
            "dots": selfies_with_dots,
            "long": long_selfies,
        }

    @pytest.mark.dependency(name="same_results")
    def test_same_results(self, datasets):
        data = [s for selfies_list in datasets.values() for s in selfies_list]
        f_old = _vectorize_and_force_computation(_split_selfies_legacy)
        f_new = _vectorize_and_force_computation(sf.split_selfies)
        assert f_old(data) == f_new(data)

    @pytest.mark.dependency(depends=["same_results"])
    @pytest.mark.parametrize("dataset_key", ["regular", "dots", "long"])
    @pytest.mark.parametrize("func", [sf.split_selfies, _split_selfies_legacy])
    def test_compare_split_methods(self, benchmark, datasets, dataset_key, func):
        benchmark.group = dataset_key
        f = _vectorize_and_force_computation(func)
        benchmark(lambda: f(datasets[dataset_key]))
        assert True
