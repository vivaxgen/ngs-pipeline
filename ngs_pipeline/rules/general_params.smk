
# necessary imports
import pathlib

# generic parameters

# get base directories
ngs_pipeline_basedir = config['NGS_PIPELINE_BASE']
ngsenv_basedir = config['NGSENV_BASEDIR']


def get_abspath(p, prefix=ngsenv_basedir):
    if p is None:
        return p

    filepath = p.as_posix() if isinstance(p, pathlib.Path) else p
    if (
        filepath.startswith("/")
        or filepath.startswith("./")
        or filepath.startswith("../")
    ):
        return p

    prefix = pathlib.Path(prefix) if isinstance(prefix, str) else prefix
    return (prefix / p).absolute().as_posix()


# basic parameters to do processing
refseq = get_abspath(config['refseq_file'], ngsenv_basedir)
refmap = get_abspath(config['refmap_file'], ngsenv_basedir)
strtable_file = get_abspath(fn if (fn := config.get('strtable_file')) else None, ngsenv_basedir)
targetregion_file = get_abspath(fn if (fn := config.get('targetregion_file')) else None, ngsenv_basedir)
ploidy = int(config.get('ploidy', 2))
min_depth = config.get('min_depth', 5)

# thread allocations
thread_allocations = config.get('thread_allocations', {})

# wildcard constraints:
wildcard_constraints:
    sample = r'[.\w-]+',
    idx = r'\d',

# for indexing, use the correct extension for bwa / bwa-mem2

bwa_bin = config.get('bwa_bin', 'bwa-mem2')
if bwa_bin == 'bwa-mem2':
    idx_extension = 'bwt.2bit.64'
elif bwa_bin == 'bwa':
    idx_extension = 'bwt'
else:
    raise RuntimeError(f'ERR: the bwa_bin {bwa_bin} is not recognized')


# for micromamba sub-environment

shell_prefix_micromamba = 'eval "$(micromamba shell hook --shell bash)" && '
shell_prefix_bowtie2 = shell_prefix_micromamba + 'micromamba activate vvg-ngspl-bowtie2 && '
shell_prefix_clair3 = shell_prefix_micromamba + 'micromamba activate vvg-ngspl-clair3 && '

# for annotation processes

gff_file = get_abspath(config["gff_file"], ngsenv_basedir)
snpEff_config_file = get_abspath(config["snpEff_config_file"], ngsenv_basedir)
snpEff_data_dir = get_abspath(config["snpEff_data_dir"], ngsenv_basedir)
snpEff_db = config["snpEff_db"]


# EOF
