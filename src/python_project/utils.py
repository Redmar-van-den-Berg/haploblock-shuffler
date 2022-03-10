def addition(a, b):
    """ Add two numbers together

    >>> addition(1, 1)
    2
    >>> addition(0, 0)
    0
    >>> addition(-1, 1)
    0
    >>> addition(-1, -1)
    -2
    """
    return a + b


def is_homozygous(call):
    """Determine if a call is homozygous"""
    return call['GT'][0] == call['GT'][1]


def is_heterozygous(call):
    """Determine if a call is heterozygous"""
    return not is_homozygous(call)


def is_compatible(call1, call2):
    """Are two calls compatible"""
    # If both are phased, the phase set must match.
    if 'PS' in call1 and 'PS' in call2:
        return call1['PS'] == call2['PS']

    # Two homozygous calls are compatible
    if is_homozygous(call1) or is_homozygous(call2):
        return True

    # Two heterozygous calls are not compatibhle
    if is_heterozygous(call1) and is_heterozygous(call2):
        return False


def are_compatible(calls, call):
    """ Determine if a call is compatible with a list of calls """
    return all((is_compatible(call, c) for c in calls))


def get_call(variant):
    """ Return the call from a variant """
    return variant.samples[0]


def group_variants(variants):
    """Group compatible variants together"""

    # Store the grouped variants by phase ID
    grouped_variants = dict()

    # Internal phase ID's for unphased variants
    internal_counter = 0

    current_group = list()
    for record in variants:
        call = get_call(record)
        if are_compatible([get_call(rec) for rec in current_group], call):
            current_group.append(record)
    grouped_variants[f"internal_{internal_counter}"] = current_group

    return grouped_variants
