import types

from python_project import utils


def test_addition():
    assert utils.addition(1, 2) == 3


def test_is_hom_ref():
    assert utils.is_homozygous({'GT': (0, 0)})


def test_is_hom_alt():
    assert utils.is_homozygous({'GT': (1, 1)})


def test_is_hom_extreme():
    assert utils.is_homozygous({'GT': (99, 99)})


def test_are_compatible_empty():
    """ Every call is compatible with an empty list """
    assert utils.are_compatible([], None)


def test_are_compatible_three_homs():
    call = {'GT': (1, 1)}
    assert utils.are_compatible([call, call], call) is True


def test_is_compatible_hom():
    """ A homozygous variant is compatible with any other variant """
    call = {'GT': (1, 1)}
    assert utils.is_compatible(call, call)


def test_is_compatible_two_het():
    """Two heterozygous calls are not compatible """
    call = {'GT': (0, 1)}
    assert utils.is_compatible(call, call) is False


def test_is_compatible_hom_het():
    """A homozygous and heterozygous call are compatible"""
    hom = {'GT': (0, 0)}
    het = {'GT': (0, 1)}
    assert utils.is_compatible(hom, het)


def test_is_compatible_phase_set():
    """Two heterozygous variants from the same phase set are compatible"""
    call = {'GT': (0, 1), 'PS': 1}
    assert utils.is_compatible(call, call)


def test_is_compatible_incompatible_phase_set():
    """Two conflicting phase sets are not compatible"""
    call1 = {'GT': (0, 1), 'PS': 1}
    call2 = {'GT': (0, 1), 'PS': 2}
    assert utils.is_compatible(call1, call2) is False


def test_is_compatible_phased_hom():
    """A call with phase set is compatible with a homozygous call"""
    hom = {'GT': (0, 0)}
    call = {'GT': (0, 1), 'PS': 1}
    assert utils.is_compatible(hom, call)


def test_is_compatible_incompatible_phased_hom():
    """A homozygous call with an incompatible phase set is not compatible"""
    hom = {'GT': (0, 0), 'PS': 0}
    call = {'GT': (0, 1), 'PS': 1}
    assert utils.is_compatible(hom, call) is False


def test_group_variants():
    hom = {'GT': (0, 0)}
    # A variant has a list of samples that contains the calls
    var = types.SimpleNamespace()
    var.samples = [hom]
    assert utils.group_variants([var]) == [[var]]


def test_group_variants_compatible():
    # Homozygous variant
    hom = {'GT': (0, 0)}
    var1 = types.SimpleNamespace(samples=[hom])

    # Heterozygous variant
    het = {'GT': (1, 0)}
    var2 = types.SimpleNamespace(samples=[het])

    assert utils.group_variants([var1, var2]) == [[var1, var2]]


def test_group_variants_incompatible():
    het = {'GT': (1, 0)}
    var = types.SimpleNamespace(samples=[het])

    assert utils.group_variants([var, var]) == [[var], [var]]


def test_generate_pattern_zero():
    assert list(utils.generate_patterns(0)) == list()


def test_generate_pattern_one():
    assert list(utils.generate_patterns(1)) == [[0]]


def test_generate_pattern_two():
    assert list(utils.generate_patterns(2)) == [[0, 0], [0, 1]]


def test_generate_pattern_four():
    expected = [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 1, 1],
            [0, 1, 0, 0],
            [0, 1, 0, 1],
            [0, 1, 1, 0],
            [0, 1, 1, 1]
    ]

    assert list(utils.generate_patterns(4)) == expected
