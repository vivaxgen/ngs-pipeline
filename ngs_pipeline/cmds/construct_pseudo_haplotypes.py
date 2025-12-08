from ngs_pipeline import cerr, cexit, arg_parser
import os 

def init_argparser():
    p = arg_parser()
    p.add_argument(
        "--variants_list",
        required=True,
        help="TSV/Bed file containing list of variants (Chrom, Pos0, Pos, Marker) for constructing haplotypes",
    )
    p.add_argument(
        "--log", default=None, help="Log file, default = None",
    )
    p.add_argument(
        "--sample",
        default="SAMPLE",
        help="Sample ID, default = SAMPLE",
    )
    p.add_argument("-o", "--outfile", required = True, help="Output filename")
    p.add_argument("infile", help="Input BAM file (collated/group by name)")
    return p

# FASTP default
WARNED_NO_PAIR_INFO = False
GOOD_QUAL = 30
BAD_QUAL = 14

def check_alignment_covers_variants(aln, aln_mate, variants):
    if aln.reference_name != aln_mate.reference_name:
        return False
    chrom = aln.reference_name
    block_start, block_end = min(aln.reference_start, aln_mate.reference_start), max(aln.reference_end, aln_mate.reference_end)
    variants_in_block = variants.query("chrom == @chrom and pos0 >= @block_start and pos < @block_end")
    if variants_in_block.shape[0] == 0:
        return False, None
    else:
        # split according to marker & pull out the whole block
        variants_in_block_grouped = []
        for marker in variants_in_block['marker'].unique():
            vars_in_marker = variants.query("marker == @marker")
            variants_in_block_grouped.append(vars_in_marker)
        return True, variants_in_block_grouped

def check_pair_info_fix(aln, aln_mate):
    if not aln.is_read1 and not aln.is_read2:
        global WARNED_NO_PAIR_INFO
        # Warn only once
        if not WARNED_NO_PAIR_INFO:
            WARNED_NO_PAIR_INFO = True
            cerr("Reads are not properly marked as read1/read2, considering forward as read1, reverse as read2")
        aln.is_read1 = aln.is_forward
        aln.is_read2 = aln.is_reverse
        aln_mate.is_read1 = aln_mate.is_forward
        aln_mate.is_read2 = aln_mate.is_reverse

    if aln.is_read1 and aln_mate.is_read2:
        return aln, aln_mate
    elif aln.is_read2 and aln_mate.is_read1:
        return aln_mate, aln
    else:
        raise ValueError("Cannot determine read1/read2 in pair")

def get_base_and_qual_at_ref_pos(aln, ref_pos):
    # This return only the alignment segment
    ref_positions = aln.get_reference_positions() 
    start, end = ref_positions[0], ref_positions[-1]
    if (ref_pos < start) or (ref_pos > end):
        return "?", -1 # not covered
    try:
        read_pos = ref_positions.index(ref_pos)
        base = aln.query_alignment_sequence[read_pos]
        qual = aln.query_alignment_qualities[read_pos]
        return base, qual
    except ValueError:
        return "-", 0 # covered, but deletion

def build_haplotype_from_variants(aln, aln_mate, covered_variants):
    try:
        r1, r2 = check_pair_info_fix(aln, aln_mate)
    except Exception as e:
        print(f"# haplotype_construction | {aln.query_name}: cannot determine read1/read2, skipping")
        return None

    haplotype = []
    marker = covered_variants['marker'].iloc[0]
    for _, var in covered_variants.iterrows():
        pos = var['pos0']
        base1, qual1 = get_base_and_qual_at_ref_pos(r1, pos)
        base2, qual2 = get_base_and_qual_at_ref_pos(r2, pos)
        
        match base1, base2:
            case "?", "?":
                haplotype.append("?")
            case b1, "?":
                haplotype.append(b1)
            case "?", b2:
                haplotype.append(b2)
            case b1, b2 if b1 == b2:
                haplotype.append(b1)
            case b1, b2 if b1 != b2:
                if qual1 >= GOOD_QUAL and qual2 <= BAD_QUAL:
                    haplotype.append(b1)
                elif qual2 >= GOOD_QUAL and qual1 <= BAD_QUAL:
                    haplotype.append(b2)
                else:
                    # prefer read1
                    haplotype.append(b1)
    print(f"# haplotype_construction | {marker}: {''.join(haplotype)} | alns: {r1.query_name}")
    return marker, "".join(haplotype)


def construct_pseudo_haplotypes(variants, alignments):
    from collections import Counter
    haplotypes = []
    tagged = 0 # alignments with tags
    total = 0  # alignments processed
    prev_aln = None

    while True:
        aln = next(alignments, None)
        print(f"Processing alignment: {total}", end="\r", flush=True, file=os.sys.stderr)
        match aln, prev_aln:
            case None, _:
                break
            case _, None:
                total += 1
                prev_aln = aln
                continue
            case aln, prev_aln if (
                aln.query_name == prev_aln.query_name and
                aln.reference_name == prev_aln.reference_name):
                total += 1
                match check_alignment_covers_variants(prev_aln, aln, variants):
                    case False, _:
                        print(f"check_pair | {prev_aln.query_name} | no interested variants")
                        prev_aln = None
                        continue
                    case True, covered_blocks:
                        print(f"check_pair | {aln.query_name} | {len(covered_blocks)} interested variant blocks")
                        for covered_variants in covered_blocks:
                            haplotype = build_haplotype_from_variants(prev_aln, aln, covered_variants)
                            if haplotype is not None:
                                haplotypes.append(haplotype)
                        tagged += 1
                        prev_aln = None
                        continue
            case _: 
                total += 1
                print(f"check_pair | {prev_aln.query_name} | mate not following or different reference")
                prev_aln = aln
                continue
    
    hap_counts = Counter(haplotypes)
    return (hap_counts, tagged, total)

def main(args):
    import pandas as pd
    import pysam

    if args.log is not None:
        log_fh = open(args.log, "w+")
        os.sys.stdout = log_fh

    bam = pysam.AlignmentFile(args.infile, "rb")
    variants = pd.read_table(args.variants_list, sep="\t", header=None, names=["chrom", "pos0", "pos", "marker"], dtype={"chrom": str, "pos0": int, "pos": int, "marker": str})
    if variants.isnull().values.any():
        cexit("Variants input file contains N/A values, please check!")
    
    haps, read_pair_with_hap, reads_processed = construct_pseudo_haplotypes(variants, bam)
    pd.DataFrame(
        [(args.sample, marker, hap, count) for (marker, hap), count in haps.items()],
        columns=["sample", "marker", "haplotype", "count"]
    ).sort_values(["marker", "count", "haplotype"]).to_csv(args.outfile, sep="\t", index=False)

    print(f"### Total reads processed | {reads_processed}")
    print(f"### Read pairs with haplotypes constructed | {read_pair_with_hap}")
    
    if args.log is not None:
        log_fh.close()
    bam.close()