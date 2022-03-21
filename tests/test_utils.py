import types
from collections import namedtuple

from haploblock_shuffler import utils


Call = namedtuple("Call", "GT")
Phased = namedtuple("Call", ["GT", "PS"])


het = Call("0/1")
hom_ref = Call("0/0")
hom = Call("1/1")
het_phased = Phased("0|1", 1)


def test_is_hom_ref():
    assert utils.is_homozygous(hom_ref)


def test_is_hom_alt():
    call = Call("1/1")
    assert utils.is_homozygous(call)


def test_is_hom_extreme():
    call = Call("99/99")
    assert utils.is_homozygous(call)


def test_is_hom_extreme_phased():
    call = Call("99|99")
    assert utils.is_homozygous(call)


def test_are_compatible_empty():
    """Every call is compatible with an empty list"""
    assert utils.are_compatible([], None)


def test_are_compatible_three_homs():
    assert utils.are_compatible([hom_ref, hom_ref], hom_ref) is True


def test_is_compatible_hom():
    """A homozygous variant is compatible with any other variant"""
    assert utils.is_compatible(hom, hom)


def test_is_compatible_two_het():
    """Two heterozygous calls are not compatible"""
    assert utils.is_compatible(het, het) is False


def test_is_compatible_hom_het():
    """A homozygous and heterozygous call are compatible"""
    assert utils.is_compatible(hom, het)


def test_is_compatible_phase_set():
    """Two heterozygous variants from the same phase set are compatible"""
    assert utils.is_compatible(het_phased, het_phased)


def test_is_compatible_incompatible_phase_set():
    """Two conflicting phase sets are not compatible"""
    call = Phased("0/1", 2)
    assert utils.is_compatible(het_phased, call) is False


def test_is_compatible_phased_hom():
    """A call with phase set is compatible with a homozygous call"""
    assert utils.is_compatible(hom, het_phased)


def test_is_compatible_incompatible_phased_hom():
    """A homozygous call with an incompatible phase set is not compatible"""
    hom = Phased("0/0", 0)
    call = Phased("0/1", 1)
    assert utils.is_compatible(hom, call) is False


def test_group_variants():
    # A variant has a list of samples that contains the calls
    var = types.SimpleNamespace()
    var.samples = [hom_ref]
    assert utils.group_variants([var]) == [[var]]


def test_group_variants_compatible():
    # Homozygous variant
    var1 = types.SimpleNamespace(samples=[hom])

    # Heterozygous variant
    var2 = types.SimpleNamespace(samples=[het])

    assert utils.group_variants([var1, var2]) == [[var1, var2]]


def test_group_variants_incompatible():
    var = types.SimpleNamespace(samples=[het])

    assert utils.group_variants([var, var]) == [[var], [var]]


def test_group_variants_phased():
    var1 = types.SimpleNamespace(samples=[het_phased])

    het2 = Phased("1/0", 2)
    var2 = types.SimpleNamespace(samples=[het2])

    assert utils.group_variants([var1, var2]) == [[var1], [var2]]


def test_group_variants_interleafed():
    """In theory, phased blocks can be interleafed"""
    var1 = types.SimpleNamespace(samples=[het_phased])

    het2 = Phased("1/0", 2)
    var2 = types.SimpleNamespace(samples=[het2])

    # var2 falls within the phase block of var1
    variants = [var1, var2, var1]

    expected = [[var1, var1], [var2]]

    assert utils.group_variants(variants) == expected


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
        [0, 1, 1, 1],
    ]

    assert list(utils.generate_patterns(4)) == expected


def test_switch_hom_ref():
    assert utils.switch(hom_ref) == hom_ref


def test_switch_het():
    assert utils.switch(het) == Call("1/0")


def test_switch_het_phased():
    assert utils.switch(het_phased) == Phased("1|0", 1)


def test_switch_call_variant():
    var = types.SimpleNamespace(samples=[het])
    switched = utils.switch_variant(var)

    # Test if we switched the calls around
    assert utils.get_call(switched) == Call("1/0")


def test_switch_variants():
    het1 = Call("1/0")
    var1 = types.SimpleNamespace(samples=[het1])

    het2 = Call("0/1")
    var2 = types.SimpleNamespace(samples=[het2])

    switched = utils.switch_variants([var1, var2])

    # Test if we got back both variants switched
    assert utils.get_call(switched[0]).GT == "0/1"
    assert utils.get_call(switched[1]).GT == "1/0"


def test_get_phase_id():
    """Get phase ID from a group of variants"""
    var = types.SimpleNamespace(samples=[het_phased])

    assert utils.get_phase_id([var, var]) == 1


def test_get_phase_none_phased():
    var = types.SimpleNamespace(samples=[het])

    assert utils.get_phase_id([var]) is None


def test_all_combinations():
    het1 = Call("1/0")
    var1 = types.SimpleNamespace(samples=[het1])

    het2 = Call("0/1")
    var2 = types.SimpleNamespace(samples=[het2])

    # pattern 00, where both are unchanged
    it = utils.all_combinations([var1, var2], 2)
    result = next(it)
    assert [utils.get_call(var[0]).GT for var in result] == ["1/0", "0/1"]

    # pattern 01, where the second variant is switched
    result = next(it)
    assert [utils.get_call(var[0]).GT for var in result] == ["1/0", "1/0"]
