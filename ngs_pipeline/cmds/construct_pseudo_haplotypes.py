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
        "--pileup", action="store_true", help="Extract count of bases only",
    )
    p.add_argument(
        "--min_qual", type=int, default=20, help="Minimum base quality to consider for haplotype construction, default = 20",
    )
    p.add_argument(
        "--min_mapq", type=int, default=30, help="Minimum mapping quality to consider for haplotype construction, default = 30",
    )
    p.add_argument(
        "--GOOD_QUAL", type=int, default=30, help="Minimum base quality to consider a base as good quality, default = 30",
    )
    p.add_argument(
        "--BAD_QUAL", type=int, default=14, help="Maximum base quality to consider a base as bad quality, default = 14",
    )
    p.add_argument(
        "--single", action="store_true", help="Reads are single-ended, not paired-end, default = False",
    )
    p.add_argument(
        "--allow_indels", action="store_true", help="Allow indels in haplotype construction, default = False (only consider SNPs)",
    )
    p.add_argument(
        "--sep", default="", help="Separator for haplotype string, default = '' (no separator)",
    )
    p.add_argument("--pileout", default=None, help="Output filename for pileup (TSV)")
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

def check_alignment_covers_variants(aln, aln_mate, variants):
    if aln_mate is None:
        aln_mate = aln
    chrom = aln.reference_name
    block_start, block_end = min(aln.reference_start, aln_mate.reference_start, key=lambda x: x if x is not None else float('inf')), \
        max(aln.reference_end, aln_mate.reference_end, key=lambda x: x if x is not None else float('-inf'))
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

def get_base_and_qual_at_ref_pos(aln, ref_pos, allow_indels=False):
    # This return only the alignment segment
    # ref_positions = aln.get_reference_positions(full_length=True) # incorrect otherwise
    # start, end = min(ref_positions, key=lambda x: x if x is not None else float('inf')), max(ref_positions, key=lambda x: x if x is not None else float('-inf'))
    aligned_pairs = aln.get_aligned_pairs(with_seq = True)
    ref_positions = [ref_pos for _, ref_pos, _ in aligned_pairs]
    if not ref_pos in ref_positions:
        return "?", -1 # not covered
    interested_pairs = [
        (read_pos, rpos, base)
        for read_pos, rpos, base in aligned_pairs
        if rpos == ref_pos
    ]
    if len(interested_pairs) == 0:
        # Distinguish unexpected no-pair cases from genuinely uncovered positions.
        cerr(f"No aligned pair found at covered reference position {ref_pos}")
        return "?", 0
    if len(interested_pairs) > 1 and not allow_indels:
        cerr(f"Multiple aligned pairs found for reference position {ref_pos}")
        return "?", 0
    elif len(interested_pairs) > 1 and allow_indels:
        bases, quals = [], []
        for read_pos, _, reference_base in interested_pairs:
            base = aln.query_sequence[read_pos]
            if base.upper() != reference_base.upper():
                print(f"Substitution: from {reference_base.upper()} to {base} at read position {ref_pos}")
            qual = aln.query_qualities[read_pos]
            bases.append(base)
            quals.append(qual)
        return bases, quals

    read_pos = interested_pairs[0][0]
    reference_base = interested_pairs[0][2]
    if read_pos is None:
        return "-", 0 # covered, but deletion
    base = aln.query_sequence[read_pos]
    if base.upper() != reference_base.upper():
        print(f"Substitution: from {reference_base.upper()} to {base} at read position {ref_pos}")
    qual = aln.query_qualities[read_pos]
    return base, qual

def build_haplotype_from_variants_single(aln, covered_variants, min_qual, allow_indels=False, sep = ""):
    haplotype = []
    marker = covered_variants['marker'].iloc[0]
    for _, var in covered_variants.iterrows():
        pos = var['pos0'] ### pos0 is 0-based, get_base_and_qual_at_ref_pos pysam use 0-based
        base, qual = get_base_and_qual_at_ref_pos(aln, pos, allow_indels)
        if allow_indels:
            if len(base) == 1:
                if base == "-":
                    haplotype.append(base)
                elif qual < min_qual:
                    haplotype.append("?")
                else:
                    haplotype.append(base)
            else:
                base_str = ""
                for b, q in zip(base, qual):
                    if b == "-":
                        base_str += b
                    elif q < min_qual:
                        base_str += "?"
                    else:
                        base_str += b
                haplotype.append(base_str)
        else:
            if qual < min_qual:
                base = "?"
            haplotype.append(base)
    print(f"# haplotype_construction_single | {marker}: {''.join(haplotype)} | aln: {aln.query_name}")
    return marker, sep.join(haplotype)

def build_haplotype_from_variants(aln, aln_mate, covered_variants, min_qual, GOOD_QUAL, BAD_QUAL, sep = ""):
    try:
        r1, r2 = check_pair_info_fix(aln, aln_mate)
    except Exception as e:
        print(f"# haplotype_construction | {aln.query_name}: cannot determine read1/read2, skipping")
        return None

    haplotype = []
    marker = covered_variants['marker'].iloc[0]
    for _, var in covered_variants.iterrows():
        pos = var['pos0'] ### pos0 is 0-based, get_base_and_qual_at_ref_pos pysam use 0-based
        base1, qual1 = get_base_and_qual_at_ref_pos(r1, pos)
        if qual1 < min_qual:
            base1 = "?"
        base2, qual2 = get_base_and_qual_at_ref_pos(r2, pos)
        if qual2 < min_qual:
            base2 = "?"

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
    return marker, sep.join(haplotype)

def deconstruct_haplotypes_alleles_to_SNPs_alleles(haplotypes, counts= None, sep = ""):
    import numpy as np
    if sep != "":
        haplotypes = [hap.split(sep) for hap in haplotypes]
    if counts is not None:
        counts = np.array(counts)
    else:
        counts = np.ones(len(haplotypes), dtype=int)
    
    observed_haplotypes = np.array([list(a) for a in haplotypes])
    assert observed_haplotypes.shape[0] == counts.shape[0]
    n_snps = observed_haplotypes.shape[1]
    alleles = []
    for snp_idx in range(n_snps):
        snp_alleles = {}
        for hap_idx in range(observed_haplotypes.shape[0]):
            allele = observed_haplotypes[hap_idx, snp_idx]
            allele_count = counts[hap_idx]
            if allele != "?":
                snp_alleles[allele] = snp_alleles.get(allele, 0) + allele_count
        alleles.append(snp_alleles)
    return alleles

def construct_pseudo_haplotypes(variants, alignments, min_mapq=30, min_qual=20, GOOD_QUAL=30, BAD_QUAL=14, pileup=False, single=False, allow_indels=False, sep=""):
    from collections import Counter
    haplotypes = []
    tagged = 0 # alignments with tags
    total = 0  # alignments processed
    prev_aln = None

    while True:
        aln = next(alignments, None)
        print(f"Processing alignment: {total}", end="\r", flush=True, file=os.sys.stderr)
        match single:
            case True:
                if aln is None:
                    break
                if aln.mapping_quality < min_mapq:
                    print(f"check_single | {aln.query_name} | low mapping quality, skipping")
                    continue
                total += 1
                match check_alignment_covers_variants(aln, None, variants):
                    case False, _:
                        print(f"check_single | {aln.query_name} | no interested variants")
                        continue
                    case True, covered_blocks:
                        print(f"check_single | {aln.query_name} | {len(covered_blocks)} interested variant blocks")
                        for covered_variants in covered_blocks:
                            haplotype = build_haplotype_from_variants_single(aln, covered_variants, min_qual, allow_indels, sep)
                            if haplotype is not None:
                                haplotypes.append(haplotype)
                        tagged += 1
                        continue
            case False:
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
                        if aln.mapping_quality < min_mapq or prev_aln.mapping_quality < min_mapq:
                            print(f"check_pair | {prev_aln.query_name} | low mapping quality, skipping")
                            prev_aln = None
                            continue
                        total += 1
                        match check_alignment_covers_variants(prev_aln, aln, variants):
                            case False, _:
                                print(f"check_pair | {prev_aln.query_name} | no interested variants")
                                prev_aln = None
                                continue
                            case True, covered_blocks:
                                print(f"check_pair | {aln.query_name} | {len(covered_blocks)} interested variant blocks")
                                for covered_variants in covered_blocks:
                                    haplotype = build_haplotype_from_variants(prev_aln, aln, covered_variants, min_qual, GOOD_QUAL, BAD_QUAL)
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
    pileup_df = variants.copy(deep =True)
    pileup_df.loc[:, ("alleles", "counts")] = "", ""
    if pileup:
        sorted_hap_counts = sorted(hap_counts.items())
        for marker in pileup_df['marker'].unique():
            haplotype_SNP_positions = pileup_df[pileup_df['marker'] == marker]
            alleles = []
            counts = []
            alleles = [haps[0][1]  for haps in sorted_hap_counts if haps[0][0] == marker]
            counts = [haps[1]  for haps in sorted_hap_counts if haps[0][0] == marker]
            if len(alleles) == 0:
                continue
            snp_alleles = deconstruct_haplotypes_alleles_to_SNPs_alleles(alleles, counts)

            for idx, row in pileup_df[pileup_df['marker'] == marker].iterrows():
                snp_idx = idx - haplotype_SNP_positions.index[0]
                snp_alleles_dict = snp_alleles[snp_idx]
                sorted_alleles = sorted(snp_alleles_dict.items(), key=lambda x: x[1], reverse=True)
                alleles_str = ",".join([allele for allele, _ in sorted_alleles]) if len(sorted_alleles) > 0 else ""
                counts_str = ",".join([str(count) for _, count in sorted_alleles]) if len(sorted_alleles) > 0 else ""
                pileup_df.at[idx, 'alleles'] = alleles_str
                pileup_df.at[idx, 'counts'] = counts_str
    else:
        pileup_df = None
    return (hap_counts, tagged, total, pileup_df)

def main(args):
    import pandas as pd
    import pysam

    if (args.allow_indels and args.sep == "") and not args.single:
        cerr("Unexpected use case")
    
    if args.log is not None:
        log_fh = open(args.log, "w+")
        os.sys.stdout = log_fh

    bam = pysam.AlignmentFile(args.infile, "rb")
    variants = pd.read_table(args.variants_list, sep="\t", header=None, names=["chrom", "pos0", "pos", "marker"], dtype={"chrom": str, "pos0": int, "pos": int, "marker": str})
    if variants.isnull().values.any():
        cexit("Variants input file contains N/A values, please check!")
    
    haps, read_pair_processed, reads_processed, SNPs = construct_pseudo_haplotypes(variants, bam, min_mapq=args.min_mapq, min_qual=args.min_qual, GOOD_QUAL=args.GOOD_QUAL, BAD_QUAL=args.BAD_QUAL, pileup=args.pileup,
        single=args.single, allow_indels=args.allow_indels, sep=args.sep)
    if args.pileup and SNPs is not None:
        if args.pileout is None:
            args.pileout = args.outfile + ".pileup.tsv"
        SNPs.to_csv(args.pileout, sep="\t", index=False)

    result_df = pd.DataFrame(
        [(args.sample, marker, hap, count) for (marker, hap), count in haps.items()],
        columns=["sample", "marker", "haplotype", "count"]
    )
    
    # fill in marker if not in result_df
    missing_markers = set(variants['marker']) - set(result_df['marker'])
    if len(missing_markers) > 0:
        marker_len = variants.query("marker in @missing_markers").groupby("marker").size().to_dict()
        missing_df = pd.DataFrame([(args.sample, marker, args.sep.join("?" * len_), 0) for marker, len_ in marker_len.items()], columns=["sample", "marker", "haplotype", "count"])
        result_df = pd.concat([result_df, missing_df], ignore_index=True)
    result_df.sort_values(["marker", "count", "haplotype"]).to_csv(args.outfile, sep="\t", index=False)
    print(f"### Total reads processed | {reads_processed}")
    print(f"### Read pairs with haplotypes constructed | {read_pair_processed}")
    
    if args.log is not None:
        log_fh.close()
    bam.close()
