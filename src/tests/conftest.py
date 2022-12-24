import pytest
import pathlib
import os


@pytest.fixture
def tests_dirpath(request):
    """get path to test directory"""
    return pathlib.Path(os.path.dirname(os.path.dirname(request.node.fspath)))

@pytest.fixture
def gtf_path_list(tests_dirpath):

    sample_suffix = '_sqanti_fltr.gtf'

    samples = ['hg002', 'hg004', 'hg005']

    gtf_list = [os.path.join(tests_dirpath, 'tests', 'data', x+sample_suffix)
                for x in samples]

    return gtf_list
