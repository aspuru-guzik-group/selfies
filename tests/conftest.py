import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--trials", action="store", default=10000,
        help="number of trails for random tests"
    )


@pytest.fixture
def trials(request):
    return int(request.config.getoption("--trials"))
