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
