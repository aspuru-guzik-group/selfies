import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--trials", action="store", default=10000,
        help="number of trails for random tests"
    )
    parser.addoption(
        "--dataset_samples", action="store", default=10000,
        help="number of samples to test from the data sets"
    )


@pytest.fixture
def trials(request):
    return int(request.config.getoption("--trials"))


@pytest.fixture
def dataset_samples(request):
    return int(request.config.getoption("--dataset_samples"))
