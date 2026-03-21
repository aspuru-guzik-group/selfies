from typing import Iterator
import pytest

import selfies as sf


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


def vectorize_and_force_computation(func):
    def wrapper(xs):
        return [list(func(x)) for x in xs]
    return wrapper


@pytest.fixture(scope="session")
def datasets():
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
def test_same_results(datasets):
    data = [s for selfies_list in datasets.values() for s in selfies_list]
    f_old = vectorize_and_force_computation(_split_selfies_legacy)
    f_new = vectorize_and_force_computation(sf.split_selfies)
    assert f_old(data) == f_new(data)


@pytest.mark.dependency(depends=["same_results"])
@pytest.mark.parametrize("dataset_key", ["regular", "dots", "long"])
@pytest.mark.parametrize("func", [sf.split_selfies, _split_selfies_legacy])
def test_compare_split_methods(benchmark, datasets, dataset_key, func):
    benchmark.group = dataset_key
    f = vectorize_and_force_computation(func)
    benchmark(lambda: f(datasets[dataset_key]))
    assert True
