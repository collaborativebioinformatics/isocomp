# pylint:disable=W0401
import pyranges as pr
import pytest

from isocomp import Coordinates
from .conftest import *


def test_window(capsys):

    # constructor
    window = Coordinates.Window('chr1', 0, 100, '+')

    assert window.chr == 'chr1'
    assert window.start == 0
    assert window.end == 100
    assert window.strand == '+'
    assert len(window) == 100

    print(window)
    captured = capsys.readouterr()
    expected = 'chr1\t0\t100\t\t1000\t+\n'
    assert captured.out == expected


def test_collapse_intervals(gtf_path_list):

    # check that all transcripts in a cluster overlap one another

    # check that there does not exist a transcript that crosses the boundary

    #

    tx_ranges = Coordinates.create_comparison_windows(gtf_path_list)

    assert 42 == 42


class TestPyrangesClustering:
    # note -- this could be moved out to a fixtures.py file and made available
    # to other tests if needed
    @pytest.fixture(autouse=True)
    def setup(self):
        """Create genomic intervals for testing pyranges concat() and 
        Cluster() methods"""
        gtf_names = ['gtf1', 'gtf2', 'gtf3']
        case_1 = [
            pr.from_dict({"Chromosome": ["chr1"], "Start": [10], "End": [20]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [30], "End": [40]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [50], "End": [60]})]

        case_2 = [
            pr.from_dict({"Chromosome": ["chr1"], "Start": [10], "End": [20]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [15], "End": [25]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [50], "End": [60]})]

        case_3 = [
            pr.from_dict({"Chromosome": ["chr1"], "Start": [10], "End": [20]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [15], "End": [25]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [18], "End": [28]})]

        case_4 = [
            pr.from_dict({"Chromosome": ["chr1"], "Start": [10], "End": [50]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [15], "End": [25]}),
            pr.from_dict({"Chromosome": ["chr1"], "Start": [30], "End": [45]})]

        # the dict isn't necessary in current tests, but i left it in b/c
        # it might be useful in the future. use the .values() method to
        # extract the list for pyragnes.concat()
        self.case_1 = dict(zip(gtf_names, case_1))
        self.case_2 = dict(zip(gtf_names, case_2))
        self.case_3 = dict(zip(gtf_names, case_3))
        self.case_4 = dict(zip(gtf_names, case_4))

    def test_no_overlap(self):
        concatenated = pr.concat(self.case_1.values())
        clustered = concatenated.cluster()
        assert clustered.Cluster.max() == 3

    def test_2_overlap(self):
        concatenated = pr.concat(self.case_2.values())
        clustered = concatenated.cluster()
        assert clustered.Cluster.max() == 2

    def test_all_overlap(self):
        concatenated = pr.concat(self.case_3.values())
        clustered = concatenated.cluster()
        assert clustered.Cluster.max() == 1

    def test_disjoint_with_common_overlapping_feature(self):
        concatenated = pr.concat(self.case_4.values())
        clustered = concatenated.cluster()
        assert clustered.Cluster.max() == 1
