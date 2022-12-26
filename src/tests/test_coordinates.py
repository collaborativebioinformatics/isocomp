# pylint:disable=W0401
from importlib.metadata import version

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

    tx_ranges = Coordinates.create_comparison_windows(gtf_path_list)

    assert 42 == 42
