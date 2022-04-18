import pytest

def pytest_addoption(parser):
    parser.addoption("--file_path", action="store", default="default file_path")
    parser.addoption("--out_path", action="store", default="default o_path")
    parser.addoption("--enrich_score", action="store", default="default enrich_score")
    parser.addoption("--copy_num", action="store", default="default copy_num")
    parser.addoption("--mer_cutoff", action="store", default="default mer_cutoff")
    parser.addoption("--merlength", action="store", default="default merlength")

@pytest.fixture
def file_path(request):
    return request.config.getoption("--file_path")

@pytest.fixture
def out_path(request):
    return request.config.getoption("--out_path")

@pytest.fixture
def enrich_score(request):
    return request.config.getoption("--enrich_score")

@pytest.fixture
def copy_num(request):
    return request.config.getoption("--copy_num")

@pytest.fixture
def mer_cutoff(request):
    return request.config.getoption("--mer_cutoff")

@pytest.fixture
def merlength(request):
    return request.config.getoption("--merlength")
