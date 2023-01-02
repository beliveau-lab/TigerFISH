import pytest

def pytest_addoption(parser):
    parser.addoption("--o_file", action="store", default="default o_file")
    parser.addoption("--r_thresh", action="store", default="default r_thresh")
    parser.addoption("--pdups_p", action="store", default="default pdups_p")
    parser.addoption("--bowtie_idx", action="store", default="default bowtie_idx")
    parser.addoption("--bt2_k_val", action="store", default="default bt2_k_val")
    parser.addoption("--max_off_target_sum", action="store", default="default max_off_target_sum")
    parser.addoption("--max_pdups_binding", action="store", default="default max_pdups_binding")

@pytest.fixture
def o_file(request):
    return request.config.getoption("--o_file")

@pytest.fixture
def r_thresh(request):
    return request.config.getoption("--r_thresh")

@pytest.fixture
def pdups_p(request):
    return request.config.getoption("--pdups_p")

@pytest.fixture
def bowtie_idx(request):
    return request.config.getoption("--bowtie_idx")

@pytest.fixture
def bt2_k_val(request):
    return request.config.getoption("--bt2_k_val")

@pytest.fixture
def max_off_target_sum(request):
    return request.config.getoption("--max_off_target_sum")

@pytest.fixture
def max_pdups_binding(request):
    return request.config.getoption("--max_pdups_binding")
