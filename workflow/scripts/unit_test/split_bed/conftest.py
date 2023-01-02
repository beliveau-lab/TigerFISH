import pytest

def pytest_addoption(parser):
    parser.addoption("--bed_file", action="store", default="default bed_file")
    parser.addoption("--chrom", action="store", default="default chrom")
    parser.addoption("--bed_out", action="store", default="default bed_out")


@pytest.fixture
def bed_file(request):
    return request.config.getoption("--bed_file")

@pytest.fixture
def chrom(request):
    return request.config.getoption("--chrom")

@pytest.fixture
def bed_out(request):
    return request.config.getoption("--bed_out")
