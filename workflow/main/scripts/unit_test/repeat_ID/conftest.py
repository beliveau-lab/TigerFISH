import pytest

def pytest_addoption(parser):
    parser.addoption("--index_file", action="store", default="default kmer_indices")
    parser.addoption("--jf_count", action="store", default="default jf_count")
    parser.addoption("--bed_file", action="store", default="default bed_file")
    parser.addoption("--chrom", action="store", default="default chrom")

@pytest.fixture
def index_file(request):
    return request.config.getoption("--index_file")

@pytest.fixture
def jf_count(request):
    return request.config.getoption("--jf_count")

@pytest.fixture
def bed_file(request):
    return request.config.getoption("--bed_file")

@pytest.fixture
def chrom(request):
    return request.config.getoption("--chrom")


