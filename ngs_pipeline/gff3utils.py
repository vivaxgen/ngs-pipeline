# gff3utils.py - ngs-pipeline library
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


def gff3_to_df(infile):
    """ return a pandas Dataframe for the 8-column GFF3 """

    import pandas as pd

    regions = {}

    # see https://m.ensembl.org/info/website/upload/gff3.html
    cols = dict(seqid=[],
                source=[],
                type=[],
                start=[],
                end=[],
                score=[],
                strand=[],
                phase=[],
    )

    with open(infile) as f_in:
        for line in f_in:
            if line.startswith('##sequence-region'):
                tokens = line.split()
                regions[tokens[1]] = int(tokens[3])
                continue

            if line.startswith('#'):
                continue

            columns = line.split('\t')

            cols['seqid'].append(columns[0])
            cols['source'].append(columns[1])
            cols['type'].append(columns[2])
            cols['start'].append(int(columns[3]))
            cols['end'].append(int(columns[4]))
            cols['score'].append(columns[5])
            cols['strand'].append(columns[6])
            cols['phase'].append(columns[7])

    df = pd.DataFrame(cols)
    # we sort by start position and longer segment
    df.sort_values(by=['seqid', 'start', 'end'],
                   ascending=[True, True, False],
                   inplace=True)

    return (regions, df)


def gff3_to_spacer(infile):

    import pandas as pd
    import dataclasses

    @dataclasses.dataclass
    class Segment:
        chrom: str
        start: int
        end: int

    merged_cols = dict(chrom=[], start=[], end=[])

    regions, df = gff3_to_df(infile)

    # we are doing a 2-pass scanning by merging segments first
    # and then find spacer based on the merged segments

    # merging

    current_segment = None

    for i, r in df.iterrows():

        chrom = r['seqid']
        start = r['start']
        end = r['end']

        if current_segment is None:
            current_segment = Segment(chrom, start, end)
            continue

        # check if we move to another chrom
        if chrom != current_segment.chrom:
            merged_cols['chrom'].append(current_segment.chrom)
            merged_cols['start'].append(current_segment.start)
            merged_cols['end'].append(current_segment.end)
            current_segment = Segment(chrom, start, end)
            continue

        # check if we are not before current segment
        if start < current_segment.start:
            raise RuntimeError('ERR: unsorted by position!')
        
        # check if we are intersect
        if start <= current_segment.end:
            if end > current_segment.end:
                # we intersect, update current segment
                current_segment.end = end
            continue

        # we are distict from current segment, update merged cols
        merged_cols['chrom'].append(current_segment.chrom)
        merged_cols['start'].append(current_segment.start)
        merged_cols['end'].append(current_segment.end)
        current_segment = Segment(chrom, start, end)
        continue

    # last segment
    merged_cols['chrom'].append(current_segment.chrom)
    merged_cols['start'].append(current_segment.start)
    merged_cols['end'].append(current_segment.end)

    merged_df = pd.DataFrame(merged_cols)

    # scanning
    spacer_cols = dict(chrom=[], start=[], end=[])

    spacer_chrom = None
    spacer_start = -1

    for i, r in merged_df.iterrows():

        chrom = r['chrom']
        start = r['start']
        end = r['end']

        if spacer_chrom is None:
            spacer_chrom = chrom
            spacer_cols['chrom'].append(spacer_chrom)
            spacer_cols['start'].append(spacer_start)
            spacer_cols['end'].append(start)
            spacer_start = end
            continue

        if chrom != spacer_chrom:
            spacer_cols['chrom'].append(spacer_chrom)
            spacer_cols['start'].append(spacer_start)
            spacer_cols['end'].append(-1)
            spacer_chrom = chrom
            spacer_start = -1
            continue

        spacer_cols['chrom'].append(spacer_chrom)
        spacer_cols['start'].append(spacer_start)
        spacer_cols['end'].append(start)
        spacer_start = end

    spacer_cols['chrom'].append(spacer_chrom)
    spacer_cols['start'].append(spacer_start)
    spacer_cols['end'].append(-1)

    spacer_df = pd.DataFrame(spacer_cols)

    # remove first and last spacer
    spacer_df.drop(
        spacer_df[(spacer_df.start < 0) | (spacer_df.end < 0)].index,
        inplace=True
    )

    spacer_df['length'] = spacer_df['end'] - spacer_df['start']
    spacer_df['midpos'] = (spacer_df['start'] + spacer_df['end'])/2

    return regions, merged_df, spacer_df


# EOF
