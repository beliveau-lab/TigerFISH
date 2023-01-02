import pytest

def pytest_addoption(parser):
    parser.addoption("--fa_file", action="store", default="default fa_file")
    parser.addoption("--chrom", action="store", default="default chrom")
    parser.addoption("--scaffold_fa", action="store", default="default scaffold_fasta out")
    parser.addoption("--jf_idx", action="store", default="default jellyfish index")
    parser.addoption("--jf_out", action="store", default="default jellyfish count output")
    parser.addoption("--k_mer_length", action="store", default="default jellyfish k-mer length")
    parser.addoption("--index_out", action="store", default="default index_out")

@pytest.fixture
def fa_file(request):
    return request.config.getoption("--fa_file")

@pytest.fixture
def chrom(request):
    return request.config.getoption("--chrom")

@pytest.fixture
def scaffold_fa(request):
    return request.config.getoption("--scaffold_fa")

@pytest.fixture
def jf_idx(request):
    return request.config.getoption("--jf_idx")

@pytest.fixture
def jf_out(request):
    return request.config.getoption("--jf_out")

@pytest.fixture
def k_mer_length(request):
    return request.config.getoption("--k_mer_length")

@pytest.fixture
def index_out(request):
    return request.config.getoption("--index_out")
