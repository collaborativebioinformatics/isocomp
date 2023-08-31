# pylint:disable=W0401
from importlib.metadata import version

from isocomp import Coordinates
from isocomp import Compare
from .conftest import *


def test_version():
    assert version('isocomp') == '0.2.0'