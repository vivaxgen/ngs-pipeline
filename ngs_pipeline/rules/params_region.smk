# parameters for defining regions

# set complete region to perform sample variant calling on all chromosome as
# single process and not per-chromosome process
# this is useful if targetregion_file is defined
complete_region = config.get('complete_region', None)

PARTIALS = None
REGIONS = config.get('regions', [])
if isinstance(REGIONS, dict):
    PARTIALS = REGIONS
    REGIONS = PARTIALS.keys()

CONTAMINANT_REGIONS = config.get('contaminant_regions', [])

wildcard_constraints:
    idx = r'\d+',
    sample = r'[.\w-]+',    # dot, alphanumerics, underscore and dash
    reg = r'[.\w-]+',       # dot, alphanumerics, underscore and dash


class RegPartition(object):

    def __init__(self, partial_regions):
        self.split = False
        self.reg_arg_name = '-L'
        self.tgt_arg_name = '-L'
        self.partitions = {}

        # sanity check
        if type(partial_regions) == dict:
            for reg, positions in partial_regions.items():
                if type(positions) != list:
                    break
                intervals = []
                start_pos = 1
                for pos in positions:
                    if pos < start_pos:
                        raise RuntimeError(
                            f'ERR: in region {reg}, {pos} < {start_pos}!!'
                        )
                    intervals.append((start_pos, pos))
                    start_pos = pos + 1
                self.partitions[reg] = intervals
            else:
                self.split = True

    def get_all_region_vcf(self, w):
        if not self.split:
            raise RuntimeError('Region is not split!')
        return list(
            [f"{destdir}/split/{w.reg}~{idx}.vcf.gz"
             for idx in range(len(self.partitions[w.reg]))
             ]
        )

    def set_arg_name(self, reg_arg_name, tgt_arg_name):
        ''' set the command line argument name for setting region
            and bed filename
        '''
        self.reg_arg_name = reg_arg_name
        self.tgt_arg_name = tgt_arg_name

    @property
    def region_vcf(self):
        if self.split:
            return temp(f"{destdir}/split/{{reg}}~{{idx}}.vcf.gz")
        return f"{destdir}/vcfs/{{reg}}.vcf.gz"

    @property
    def notation(self):
        if self.split:
            return '{reg}~{idx}'
        return '{reg}'

    def get_interval(self, w):
        if self.split:
            # we will have w.reg and w.idx
            if not hasattr(w, 'idx'):
                raise RuntimeError('wildcards do not have idx variable')
            idx = int(w.idx)
            start, end = self.partitions[w.reg][idx]
            return f'{self.reg_arg_name} {w.reg}:{start}-{end}'

        if interval_dir:
            return f'{self.tgt_arg_name} {interval_dir}/{w.reg}.bed'

        if interval_file:
            return f'{self.tgt_arg_name} {interval_file}'

        if complete_region:
            return " ".join(f'{self.reg_arg_name} {reg}' for reg in REGIONS)

        return f'{self.reg_arg_name} {w.reg}'


# EOF
