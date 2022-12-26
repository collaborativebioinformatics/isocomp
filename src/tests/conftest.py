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


@pytest.fixture
def clustered_gtf(tests_dirpath):

    tests_dirpath = os.path.join(tests_dirpath, 'tests', 'data')

    return os.path.join(tests_dirpath, 'clustered_regions.gtf')


@pytest.fixture
def fasta_dict(tests_dirpath):

    tests_dirpath = os.path.join(tests_dirpath, 'tests', 'data')

    d = dict(zip(['hg002_sqanti_fltr',
                  'hg004_sqanti_fltr',
                  'hg005_sqanti_fltr'],
                 [os.path.join(tests_dirpath, 'hg002_sqanti_fltr.fasta'),
                  os.path.join(tests_dirpath, 'hg004_sqanti_fltr.fasta'),
                  os.path.join(tests_dirpath, 'hg005_sqanti_fltr.fasta')]))

    return d
