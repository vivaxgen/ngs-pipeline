from ngs_pipeline import cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser()
    p.add_argument(
        "--phased-vcf", required=True, help="phased VCF to use for recode haplotypes"
    )
    p.add_argument("-o", "--outfile")
    p.add_argument("infile")

    return p


def recode_haplotypes(args):

    import cyvcf2

    # import pandas as pd

    # parse VCF file
    vcf = cyvcf2.VCF(args.phased_vcf)

    variant_positions = list(vcf())

    # read tsv file
    with open(args.infile) as f_in:
        haplotypes = [s.split() for s in f_in]

    recoded_haplotypes = []
    for haplotype, count, ratio, sample in haplotypes:

        pos_alleles = haplotype.split("|")
        recoded_pos_alleles = {}
        recoded_positions = set()

        for pos_allele, v in zip(pos_alleles, variant_positions):
            pos, allele = pos_allele.split(":")

            if int(pos) != v.POS:
                raise ValueError(f"No matching: {pos=} and {v.POS=}")

            if v.POS in recoded_positions:
                continue

            bcsq_info = v.INFO.get("BCSQ", None)
            if not bcsq_info:
                raise ValueError(f"No BCSQ in {v.POS=}")

            if bcsq_info.startswith("@"):
                continue

            token = bcsq_info.split("|")
            codons = token[-2].split(">")
            positions = [int(p[:-3]) for p in token[-1].split("+")]

            allele_idx = 0 if allele == v.REF else 1

            cdn = codons[allele_idx]
            cdn_pos, aa = int(cdn[:-1]), cdn[-1]
            recoded_pos_alleles[cdn_pos] = aa
            recoded_positions |= set(positions)

        aa_alleles = [f"{p}:{a}" for (p, a) in sorted(recoded_pos_alleles.items())]
        recoded_haplotypes.append(("|".join(aa_alleles), count, ratio, sample))

    with open(args.outfile, "w") as f_out:
        for items in recoded_haplotypes:
            f_out.write("\t".join(items) + "\n")


def main(args):
    recode_haplotypes(args)


# EOF
