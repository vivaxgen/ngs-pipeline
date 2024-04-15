__copyright__ = '''
prepare_samples.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2022-2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under MIT license.
Please read the README.txt of this software.
'''

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
from ngs_pipeline import cerr, cexit, arg_parser, check_NGSENV_BASEDIR


def init_argparser():
    p = arg_parser(desc='prepare directory structure for sample processing')
    p.add_argument('--resampling', type=int, default=-1,
                   help='randomly take a number of samples')
    p.add_argument('-o', '--outdir', default='analysis',
                   help='the name of the output directory, default to "analysis"')
    p.add_argument('-i', '--infile', default='-',
                   help='the name of file containing sample and read manifest, '
                        'default to stdin')
    p.add_argument('-f', '--force', default=False, action='store_true',
                   help='force the preparation even if the output directory is not '
                        'under current pipeline base directory')
    p.add_argument('indir',
                   help='the name of the directory containing reads')
    return p


def prepare_samples(args):

    # this will create the directory structures as following:
    # outdir/[SAMPLE]/reads/
    #                       raw-0_R1.fastq.gz
    #                       raw-0_R2.fastq.gz
    # each .fastq.gz is a symbolic link to real file

    # import heavy modules here if required
    import pathlib
    import pandas as pd

    # read manifest file
    in_stream = sys.stdin if args.infile == '-' else args.infile
    manifest_df = pd.read_table(in_stream, sep=None, engine='python')

    # get absolute path for indir and outdir, without resolving symlinks
    indir = pathlib.Path(args.indir).absolute()
    outdir = pathlib.Path(args.outdir).absolute()

    if not args.force:
        NGSENV_BASEDIR = check_NGSENV_BASEDIR()

        # check whether outdir is under NGSENV_BASEDIR
        if not outdir.resolve().is_relative_to(NGSENV_BASEDIR):
            cexit(f'[ERROR: {outdir} is not relative to {NGSENV_BASEDIR}]')

    # check and search for indir & outdir common parent directory
    common_parent = None
    for idx in range(len(outdir.parents)):
        if indir.is_relative_to(outdir.parents[idx]):
            common_parent = outdir.parents[idx]
            break

    if not common_parent:
        cexit(f'Directory {args.indir} and {args.outdir} do not have common '
              f'parent directory')

    # create relative link reference
    rel_path = indir.relative_to(common_parent)
    link_path = ['..'] * (idx + 3) + [rel_path]
    source_dir = pathlib.Path(*link_path)
    dest_dir = outdir.relative_to(common_parent)

    if args.resampling > 0:
        manifest_df = manifest_df.sample(n=args.resampling)

    # iterating over manifest file and check ooccurence of each fastq file
    counter = 0
    samples = []
    for (idx, r) in manifest_df.iterrows():
        sample = r['SAMPLE']
        reads = r['FASTQ']

        if type(reads) != str or type(sample) != str:
            cexit(f'ERROR: invalid values in the row with SAMPLE: {sample} '
                  f'and FASTQ: {reads}',
                  err_code=11)

        sample = sample.strip()
        reads = reads.strip()

        if not reads or not sample:
            cexit(f'ERROR: missing values in the row with SAMPLE: {sample} '
                  f'and FASTQ: {reads}',
                  err_code=12)

        if sample.startswith('#'):
            continue

        fastq_list = []
        # split reads for multiple runs
        for fastq_pair in reads.split(';'):

            path_pair = []
            for fastq_file in fastq_pair.split(','):
                fastq_path = indir / fastq_file
                if not fastq_path.is_file():
                    cexit(f'ERROR: path {fastq_path} does not exist. Plase check '
                          f'manifest file line {idx+1}')
                path_pair.append(fastq_file)

            fastq_list.append(path_pair)

        samples.append((sample, fastq_list))
        counter += 1

    # preparing directory structure

    outdir.mkdir(exist_ok=True)

    # sanity check for duplicate sample (directory) name
    duplicated = []
    for (sample, fastq_list) in samples:
        sample_path = dest_dir / sample
        if sample_path.exists():
            duplicated.append(sample)
    if any(duplicated):
        cerr('ERROR: the directories for the following samples already exist:')
        for c in duplicated:
            cerr(f'  {c}')
        cexit('Please either remove the directories or remove the sample from manifest file!')

    # for each samples, create a directory reads
    for (sample, fastq_list) in samples:

        cerr(f'Preparing for sample [{sample}]')
        sample_path = dest_dir / sample / 'reads'
        sample_path.mkdir(parents=True)

        for idx, fastq_pair in enumerate(fastq_list):
            for no, fastq_file in enumerate(fastq_pair, 1 if len(fastq_pair) > 1 else 0):
                dest = sample_path / f'raw-{idx}_R{no}.fastq.gz'
                dest.symlink_to(source_dir / fastq_file)

    cerr(f'Finished preparing {len(samples)} sample(s).')


def main(args):
    prepare_samples(args)

# EOF
