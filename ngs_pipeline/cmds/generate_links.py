__copyright__ = '''
generate_links.py - ngs-pipeline command line
[https://github.com/vivaxgen/ngs-pipeline]

(c) 2023 Hidayat Trimarsanto <trimarsanto@gmail.com>

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
    p = arg_parser(desc='generate symbolic links for each of input files')
    p.add_argument('-o', '--outdir', default='analysis',
                   help='the name of the output directory [infiles]')
    p.add_argument('-u', '--underline', type=int, default=0,
                   help='the number of underline for filename splitting, counted in reverse')
    p.add_argument('--use-absolute-link', default=False, action='store_true',
                   help='Use absolute link instead of relative link')
    p.add_argument('infiles', nargs='+',
                   help='input files in compressed FASTQ format (.fastq.gz)')
    return p


def generate_links(args):

    # import heavy modules here if required
    import pathlib
    from ngs_pipeline import fileutils

    # generate read files dictionary
    read_files = fileutils.ReadFileDict(args.infiles, args.underline)

    if any(read_files.err_files):
        cexit('ERROR: invalid input files: \n' +
            '\n'.join(f'  {errmsg}' for errmsg in read_files.err_files))

    # for each sample, process file
    counter = 0
    for sample in read_files.samples():
        # check if the current sample has more than 1 read file (or paired files)
        items = read_files[sample]
        if len(items) > 1:
            # accommodate by addding ~idx (tilde), eg: sample-01~0_R1.fastq.gz
            for idx, item in enumerate(items):
                fileutils.make_sample_symlink(sample, item, args.outdir, idx,
                                              use_absolute=args.use_absolute_link)
        else:
            fileutils.make_sample_symlink(sample, items[0], args.outdir,
                                          use_absolute=args.use)
        
        cerr(f'Generated symlink for sample {sample}')
        counter += 1

    cerr(f'Generated symlink for {counter} sample(s).')


def main(args):
    generate_links(args)

# EOF
