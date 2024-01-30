
# utilities for dealing with references

import yaml

__template__ = {
    'profile_code': '',
    'refseq_file': '',
    'refmap_file': '',
    'knownsite_file': '',
    'target_regions': '',
    'target_variants': '',
    'variant_info': ''
}


def get_ngs_profile(code, yamlfile):
    d = yaml.load(yamlfile)
    return d[code]

