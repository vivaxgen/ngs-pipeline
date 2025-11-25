#!/usr/bin/env ngs-pl


from ngs_pipeline import cerr, cexit, arg_parser


def init_argparser():

    p = arg_parser()
    p.add_argument(
        "--min-match-length", type=int, default=-1, help="minimum length of matches"
    )
    p.add_argument(
        "--min-qual", type=float, default=15.0, help="minimum value of QUAL in VCF"
    )
    p.add_argument(
        "--min-depth-ratio",
        type=float,
        default=0.5,
        help="minimum ratio of depth against median of total depths of all variants",
    )
    p.add_argument(
        "--vcf-file",
        required=True,
        help="VCF file containing list of variants for constructing haplotypes",
    )
    p.add_argument(
        "--cds-only",
        action="store_true",
        default=False,
        help="only use CDS variants for constructing haplotypes",
    )
    p.add_argument(
        "--sample",
        required=True,
    )
    p.add_argument("--outvcf", default="", help="VCF output file with phased genotype")
    p.add_argument("--outlog", default="", help="log file in YAML format")
    p.add_argument("-o", "--outfile", help="Output filename")
    p.add_argument("infile", help="Input BAM file")
    return p


def build_tag_dict(pair):
    d = {}
    for p in pair:
        k = p[1]
        if k:
            d[k] = p
    return d


def get_tag_at_ref_pos(pos, pairs):
    # pairs = [(3919, 404938, 'A', <CIGAR_OPS.CEQUAL: 7>), ...]

    import bisect

    return bisect.bisect_left(pairs, pos, key=lambda x: x[1])


def write_new_vcf(fname, vcf, variants, haplotypes):

    import cyvcf2

    # transposed haplotypes
    haplotypes_T = [list(r) for r in zip(*haplotypes)]

    w = cyvcf2.Writer(fname, vcf)

    for alleles, variant in zip(haplotypes_T, variants):
        allele_list = [variant.REF] + variant.ALT
        variant.genotypes[0] = [allele_list.index(a[-1]) for a in alleles] + [True]
        variant.genotypes = variant.genotypes
        w.write_record(variant)

    w.close()


def generate_haplotype_counter(bam, variant_positions, args):

    from collections import Counter

    alignments = bam.fetch()
    tagged = 0
    total = 0
    haplotypes = []

    for aln in alignments:
        try:
            total += 1
            paired_positions = aln.get_aligned_pairs(
                matches_only=True, with_seq=True, with_cigar=True
            )
            if (args.min_match_length > 0) and (
                len(paired_positions) < args.min_match_length
            ):
                continue
            tagged += 1

        except ValueError:
            continue

        except TypeError:
            continue

        hap = []

        for varpos in variant_positions:
            if aln.reference_name != varpos[0]:
                cerr("alignment is not part of reference")
                break

            # paired_positions use 0-based offset, while varpos[1] is 1-based
            # offset, hence decrement of varpos[1]
            pos = get_tag_at_ref_pos(varpos[1] - 1, paired_positions)

            if (pos != 0) and (pos < len(paired_positions)):
                paired_pos = paired_positions[pos]
                # paired_pos is:
                # (query_pos, ref_pos, ref_allele [lower case if substition], ops)
                # (3919, 404938, 'A', <CIGAR_OPS.CEQUAL: 7>), ...]

                base = aln.query_sequence[paired_pos[0]]

                # position checking
                if paired_pos[1] + 1 != varpos[1]:
                    cerr(f"possible deletion at {varpos[1]=}")
                    break

                # consistency checking if ref_allele equal to VCF REF
                if paired_pos[2].upper() != varpos[2]:
                    base = "^"

                if base != varpos[2] and base not in varpos[3]:
                    base = "?"

                hap.append(f"{paired_pos[1] + 1}:{base}")
            else:
                # the current read does not cover all target variants or
                # does not have all pairing base
                break

        else:
            # import IPython; IPython.embed()
            haplotypes.append(tuple(hap))

    hap_counts = Counter(haplotypes)

    return (hap_counts, tagged, total)


def construct_haplotypes(args):

    import sys
    import statistics as stats
    from collections import Counter

    import pysam
    import cyvcf2

    # read VCF file
    vcf = cyvcf2.VCF(args.vcf_file, samples=[args.sample])

    variant_positions = []
    variants = list(vcf())
    for v in variants:
        if v.QUAL < args.min_qual:
            # remove low QUAL variants
            continue
        if len(v.REF) > 1 or any([len(a) > 1 for a in v.ALT]):
            # remove non single nucleotide variants
            continue
        if args.cds_only:
            bcsq = v.INFO.get("BCSQ", None)
            if not bcsq:
                cexit(
                    "ERR: VCF file does not have BCSQ tag in INFO field; "
                    "run bcftools csq first on the VCF file"
                )
            if ";" in bcsq:
                raise NotImplementedError(
                    "multiple field in BCSQ currently not handled"
                )
            tokens = bcsq.split("|")
            if tokens[0] not in ["missense", "synonymous"]:
                continue

        variant_positions.append((v.CHROM, v.POS, v.REF, v.ALT, v.gt_depths[0], v))

    # filtering variant for total depths
    if any(variant_positions):
        mean_depth = (
            stats.mean([t[4] for t in variant_positions]) * args.min_depth_ratio
        )
        variant_positions = [v for v in variant_positions if v[4] > mean_depth]

    cerr(f"Haplotyping for {len(variant_positions)} variant(s)")

    bam = pysam.AlignmentFile(args.infile, "rb")

    hap_counts, tagged_reads, total_reads = generate_haplotype_counter(
        bam, variant_positions, args
    )

    if args.outfile:
        f_out = open(args.outfile, "w")
    else:
        f_out = sys.stdout

    reported_hap_total = 0
    reported_haps = []
    for h, count in hap_counts.items():
        proportion = count / hap_counts.total()
        if proportion < 0.15 or count < 5:
            continue

        reported_hap_total += count
        reported_haps.append(h)
        f_out.write(f"{'|'.join(h)}\t{proportion:5.3f}\t{count}\t{args.sample}\n")

    # import IPython; IPython.embed()
    if args.outlog:
        import json

        d = dict(
            sample=args.sample,
            total_reads=total_reads,
            tagged_reads=tagged_reads,
            haplotyped_reads=hap_counts.total(),
            haplotyped_reads_ratio=hap_counts.total() / tagged_reads,
            haplotyped_reported=reported_hap_total,
            haplotyped_reported_ratio=reported_hap_total / hap_counts.total(),
        )
        with open(args.outlog, "w") as json_fout:
            json.dump(d, json_fout)

    cerr(f"total reads: {total_reads}")
    cerr(f"tagged reads: {tagged_reads}")
    cerr(
        f"haplotyped reads: {hap_counts.total()} ({hap_counts.total() / tagged_reads :5.3f})"
    )
    cerr(
        f"haplotyped reported: {reported_hap_total} ({reported_hap_total / hap_counts.total() :5.3f})"
    )

    if args.outvcf:
        write_new_vcf(
            args.outvcf, vcf, [v[-1] for v in variant_positions], reported_haps
        )


def main(args):
    construct_haplotypes(args)


# EOF
