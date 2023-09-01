import pytest
import os

@pytest.fixture
def tests_dirpath(request):
    """Get the path to the test directory."""
    return os.path.dirname(os.path.dirname(request.node.fspath))

@pytest.fixture
def gtf_path_list(tests_dirpath):
    sample_suffix = '_sqanti_fltr.gtf'
    samples = ['hg002', 'hg004', 'hg005']
    gtf_list = [os.path.join(tests_dirpath, 'tests', 'data', f"{x}{sample_suffix}") for x in samples]
    return gtf_list

@pytest.fixture
def clustered_gtf(tests_dirpath):
    return os.path.join(tests_dirpath, 'tests', 'data', 'clustered_regions.gtf')

@pytest.fixture
def fasta_dict(tests_dirpath):
    tests_dir = os.path.join(tests_dirpath, 'tests', 'data')
    sample_names = ['hg002_sqanti_fltr', 'hg004_sqanti_fltr', 'hg005_sqanti_fltr']
    d = {name: os.path.join(tests_dir, f"{name}.fasta") for name in sample_names}
    return d
