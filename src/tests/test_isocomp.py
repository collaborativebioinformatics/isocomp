# pylint:disable=W0401
import logging
from importlib.metadata import version

from isocomp import Coordinates
from .conftest import *


def test_version():
    assert version('isocomp') == '0.1.0'


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
    assert captured.out == 'chr1\t0\t100\t+\n'



def test_collapse_intervals(gtf_path_list):

    tx_ranges = Coordinates.create_comparison_windows(gtf_path_list)

    assert 2 == 2
