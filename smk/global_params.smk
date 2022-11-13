
import pathlib

sample = pathlib.Path.cwd().name

platform = 'ILLUMINA' if config.get('instrument', '').lower() in ['miseq', 'nextseq', 'novaseq'] else config['platform']


