from python_project import utils


def test_addition():
    assert utils.addition(1, 2) == 3


def test_is_hom_ref():
    assert utils.is_homozygous({'GT': (0, 0)})


def test_is_hom_alt():
    assert utils.is_homozygous({'GT': (1, 1)})


def test_is_hom_extreme():
    assert utils.is_homozygous({'GT': (99, 99)})


def test_compatible_empty():
    """ Every call is compatible with an empty list """
    assert utils.compatible([], None)


def test_compatible_hom():
    """ A homozygous variant is compatible with any list of variants """
    assert utils.compatible([1], {'GT': (1, 1)})


def test_compatible_two_het():
    """Two heterozygous calls are not compatible """
    call = {'GT': (0, 1)}
    assert utils.compatible([call], call) is False


def test_compatible_hom_het():
    """A homozygous and heterozygous call are compatible"""
    hom = {'GT': (0, 0)}
    het = {'GT': (0, 1)}
    assert utils.compatible([hom], het)


def test_compatible_phase_set():
    """Two heterozygous variants from the same phase set are compatible"""
    call = {'GT': (0, 1), 'PS': 1}
    assert utils.compatible([call], call)


def test_incompatible_phase_set():
    """Two conflicting phase sets are not compatible"""
    call1 = {'GT': (0, 1), 'PS': 1}
    call2 = {'GT': (0, 1), 'PS': 2}
    assert utils.compatible([call1], call2) is False
