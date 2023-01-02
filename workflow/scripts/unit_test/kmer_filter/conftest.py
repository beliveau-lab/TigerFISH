import pytest

def pytest_addoption(parser):
    parser.addoption("--probe_file", action="store", default="default bed_file")
    parser.addoption("--fasta_file", action="store", default="default fasta_file")
    parser.addoption("--jf_file", action="store", default="default jf_file")
    parser.addoption("--o_path", action="store", default="default out_file")

@pytest.fixture
def probe_file(request):
    return request.config.getoption("--probe_file")

@pytest.fixture
def fasta_file(request):
    return request.config.getoption("--fasta_file")

@pytest.fixture
def jf_file(request):
    return request.config.getoption("--jf_file")

@pytest.fixture
def o_path(request):
    return request.config.getoption("--o_path")
