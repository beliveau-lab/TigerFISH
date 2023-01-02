import pytest

def pytest_addoption(parser):
    parser.addoption("--bed", action="store", default="default bed")
    parser.addoption("--region_fa", action="store", default="default region_fa")
    parser.addoption("--genome_fa", action="store", default="default genome_fa")
    parser.addoption("--test_regions", action="store", default="default test_regions")
    parser.addoption("--probe_out", action="store", default="default probe_out")
    parser.addoption("--name", action="store", default="default name")


@pytest.fixture
def bed(request):
    return request.config.getoption("--bed")

@pytest.fixture
def region_fa(request):
    return request.config.getoption("--region_fa")

@pytest.fixture
def genome_fa(request):
    return request.config.getoption("--genome_fa")

@pytest.fixture
def name(request):
    return request.config.getoption("--name")

@pytest.fixture
def test_regions(request):
    return request.config.getoption("--test_regions")

@pytest.fixture
def probe_out(request):
    return request.config.getoption("--probe_out")


