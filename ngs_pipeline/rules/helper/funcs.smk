
from ngs_pipeline import cerr

def _expand_sp(w, pattern):
    """ expand the given pattern with the sample prefix, sample name and/or index if necessary. """

    if "<sp>" in pattern:
        # replace <sp> with the sample prefix
        pattern = pattern.replace("<sp>", SP)

        if "{idx}" in pattern:
            # if the pattern contains {idx}, we need to expand it with the list of indexes for the given sample
            idxs = get_indexes(w)
            if "{sample}" in pattern:
                if "{pfx}" in pattern:
                    # if the pattern contains {sample}, we need to expand it with the sample name as well
                    cerr(f"expanding pattern: {pattern} with pfx: {w.pfx}, sample: {get_sample(w)} and indexes: {idxs}")
                    return expand(pattern, pfx=w.pfx, sample=get_sample(w), idx=idxs)

                cerr(f"expanding pattern: {pattern} with sample: {get_sample(w)} and indexes: {idxs}")
                return expand(pattern, sample=get_sample(w), idx=idxs)
                
            # just expand the pattern with the list of indexes for the given sample
            cerr(f"expanding pattern: {pattern} with indexes: {idxs}")
            return expand(pattern, idx=idxs)
        # if the pattern does not contain {idx}, we just return the pattern with the sample prefix and sample name
        return [pattern.format(pfx=w.pfx, sample=get_sample(w))]

    if "{idx}" in pattern:
        # if the pattern contains {idx}, we need to expand it with the list of indexes for the given sample
        idxs = get_indexes(w)
        return expand(pattern, idx=idxs)
    
    # if the pattern does not contain <sp> or {idx}, we just return the pattern as is
    return pattern


def expand_sp(pattern):
    """
    This function is used to expand the given pattern with the sample prefix, sample name and/or index if necessary.
    DO NOT use this function in the output directive (similar like expand()).
    """
    return lambda w: _expand_sp(w, pattern)


def replace_sp(pattern):
    """
    This function is used to replace <sp> to actual path.
    This function can be used in the output directive.

    """
    return pattern.replace("<sp>", SP)


def temp_unless(path, condition):
    """set the path as temporary unless the condition is True. """
    return path if condition else temp(path)


def temp_unless_config(path, config_key):
    """ set the path as temporary unless the config key is True. """
    return temp_unless(path, config.get(config_key, False))


# EOF
