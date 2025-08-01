
import pathlib

include: "general_params.smk"

# parameters to do processing

strtable_file = get_abspath(config.get('strtable_file', 'STRTABLE-FILE-NOT-DEFINED'), ngsenv_basedir)
knownsites_file = get_abspath(config.get('knownsites_file', ''), ngsenv_basedir)
knownvariants_dir = get_abspath(config.get('knownvariants_dir', ''), ngsenv_basedir)
targetregion_file = get_abspath(fn if (fn := config.get('targetregion_file')) else None, ngsenv_basedir)

deduplicate = config.get('deduplicate', True)

# parameters to save storage space
keep_paired_bam = config.get('keep_paired_bam', False)
keep_proper_bam = config.get('keep_proper_bam', False)
keep_final_bam = config.get('keep_final_bam', False)
keep_filtered_bam = config.get('keep_filtered_bam', False) or (keep_final_bam if not deduplicate else False)
keep_deduplicated_bam = config.get('keep_deduplicated_bam', False) or keep_final_bam
keep_recalibrated_bam = config.get('keep_recalibrated_bam', False)

# parameters for read trimming
minlen = int(config['minlen']) if 'minlen' in config else int(config['read_length'] / 3)
maxlen = int(config['maxlen']) if 'maxlen' in config else 0
min_read_quality = int(config.get('min_read_qual', 15))
min_avg_quality = int(config.get('min_avg_qual', 0))
correction = config.get('correction', False)
instrument = config.get('instrument', None)
read_filters = config.get('read_filters', '')

# parameter for joint variant calling
variant_file = config.get('variant_file', None)
interval_file = config.get('interval_file', None)
interval_dir = config.get('interval_dir', None)

# generic parameters
java_opts = config.get('java_opts', '')

# check available read files

ngs_platform = 'ILLUMINA' if config.get('instrument', '').lower() in ['miseq', 'nextseq', 'novaseq'] else config['platform']

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
