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


def compatible(calls, call):
    """ Determine if a call is compatible with a list of calls """
    # Any variant is compatible with an empty list of cals
    if not calls:
        return True

    # A homozygous call is compatible with any list of calls
    if is_homozygous(call):
        return True

    # Check each call in turn, does the phase set match for any call
    for c in calls:
        if 'PS' in call:
            return call['PS'] == c.get('PS')

    # If call is heterozygous, all (other) existing calls must be homozygous
    return all((is_homozygous(c) for c in calls))
