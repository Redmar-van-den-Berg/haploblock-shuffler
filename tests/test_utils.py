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
