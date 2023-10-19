#!/usr/bin/env python3

from ngsutils import cerr, arg_parser

from pathlib import Path
from PIL import Image


def init_argparser():

    p = arg_parser('gather depth plots')
    p.add_argument('-n', type=int, default=1,
                   help='number of sample(s) per page')
    p.add_argument('-p', '--plotfile', default='logs/depths.png',
                   help='name of plot file [logs/depth.png]')
    p.add_argument('-o', '--outfile', default='outplot.pdf')
    p.add_argument('indirs', nargs='+')

    return p


def batch(seq, size):
    return [
        seq[i:i + size]
        for i in range(0, len(seq), size)
    ]


def gather_depth_plots(args):

    plotfiles = []

    for indir in args.indirs:
        for plotfile in Path(indir).glob('*/' + args.plotfile):
            plotfiles.append(plotfile)

    cerr(f'[Gathering {len(plotfiles)} depth plots]')

    pages = []
    counter = 0
    for file_list in batch(plotfiles, args.n):

        imgs = [Image.open(infile) for infile in file_list]
        imgs = [img.resize((img.width//2, img.height//2)) for img in imgs]

        img_height = sum([img.height for img in imgs])
        img_width = max([img.width for img in imgs])

        dest = Image.new('RGB', (img_width, img_height))

        y_offset = 0
        for img in imgs:
            dest.paste(img, (0, y_offset))
            y_offset += img.height

        pages.append(dest)
        counter += len(imgs)

    pages[0].save(args.outfile, format='PDF', save_all=True, append_images=pages[1:])
    print('Processing depthplot for %d samples' % counter)


def main(args):
    gather_depth_plots(args)

# EOF
