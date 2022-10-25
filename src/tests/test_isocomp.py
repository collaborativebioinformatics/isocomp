from importlib.metadata import version 

def test_version():
    assert version('isocomp') == '0.1.0'