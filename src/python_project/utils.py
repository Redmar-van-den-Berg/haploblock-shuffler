import copy


def is_homozygous(call):
    """Determine if a call is homozygous"""
    return call["GT"][0] == call["GT"][1]


def is_heterozygous(call):
    """Determine if a call is heterozygous"""
    return not is_homozygous(call)


def is_compatible(call1, call2):
    """Are two calls compatible"""
    # If both are phased, the phase set must match.
    if "PS" in call1 and "PS" in call2:
        return call1["PS"] == call2["PS"]

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


def get_phase_id(variants):
    for var in variants:
        call = get_call(var)
        if 'PS' in call:
            return call['PS']


def add_group(grouped_variants, current_group, phased):
    """Add current_group to the grouped variants

    This can be tricky when there are phased variants
    """
    phase_id = get_phase_id(current_group)
    if not phase_id:
        grouped_variants.append(current_group)
    else:
        # If we have seen this phase ID before
        if phase_id in phased:
            index = phased[phase_id]
            grouped_variants[index] += current_group
        else:
            index = len(grouped_variants)
            grouped_variants.append(current_group)
            phased[phase_id] = index


def group_variants(variants):
    """Group compatible variants together"""

    # Store the grouped variants in a list of list of variants
    grouped_variants = list()

    # Store the index of phased groups, since they could be interleafed, so we
    # have to add more variants to them lates
    phased = dict()

    current_group = list()
    for record in variants:
        call = get_call(record)
        if are_compatible([get_call(rec) for rec in current_group], call):
            current_group.append(record)
        else:
            add_group(grouped_variants, current_group, phased)
            current_group = [record]

    # If we were still working on a group of variants when we got to the last
    # one
    if current_group:
        add_group(grouped_variants, current_group, phased)

    return grouped_variants


def generate_patterns(count):
    """Generate patterns for switching variants around"""
    if count < 1:
        return list()

    for i in range(2 ** (count - 1)):
        yield [int(x) for x in format(i, "b").zfill(count)]


def switch(call):
    """Switch the genotype calls around"""
    new = call.copy()
    new["GT"] = new["GT"][::-1]
    return new


def switch_variant(var):
    """ Switch the calls for a variant around """
    newvar = copy.deepcopy(var)
    call = get_call(newvar)
    newvar.samples = [switch(call)]
    return newvar


def switch_variants(variants):
    """Switch the calls for all variants"""
    return [switch_variant(var) for var in variants]


def all_combinations(variants):
    """Yield all possible combinations of variants"""
    grouped = group_variants(variants)

    for pattern in generate_patterns(len(grouped)):
        yield [
            switch_variants(group) if p else group for group, p in zip(grouped, pattern)
        ]
