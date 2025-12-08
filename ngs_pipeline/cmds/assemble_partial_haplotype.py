from ngs_pipeline import cerr, cexit, arg_parser
import os


def init_argparser():
    p = arg_parser()
    p.add_argument(
        "--log", default=None, help="Log file, default = None",
    )
    p.add_argument("-o", "--outfile", required = True, help="Output filename")
    p.add_argument("-m", "--max-consider", type=int, default=15,
        help="Maximum number of seeds, default = 15")
    p.add_argument("-p", "--prioritise", default="count",
        choices=["count", "completeness"],
        help="Prioritise haplotypes by count or completeness (number of missing bases)")
    p.add_argument("infile", help="haplotype tsv for a single sample")
    return p

def is_compatible(hap1, hap2) -> bool:
    return all((base1 == base2 or base1 == "?" or base2 == "?") for base1, base2 in zip(hap1, hap2))

def scan_for_seed(haplotypes: list[str]) -> list[str]:
    seeds = []
    for hap in haplotypes:
        if set(list(hap)) == {"?"}:
            continue
        elif hap.count("?") == len(hap) - 1:
            continue
        elif seeds == []:
            seeds.append(hap)
        else:
            has_compatible = False
            for seed in seeds:
                # if incompatible with all existing seeds, add as new seed
                if is_compatible(hap, seed):
                    has_compatible = True
                    break
            if not has_compatible:
                seeds.append(hap)
    return seeds

def assign_to_groups(hap_df, seeds) -> list[dict]:
    # seeds is predetermined, all haplotypes should be able to be assigned to one of the seeds,
    # if seed has been filtered out, sequence can be ignored
    groups = {s: {"h": [], "t":[]} for s in seeds}

    for _, row in hap_df.iterrows():
        haplotype = row["haplotype"]
        count = row["count"]
        # remove haplotype with no interestde SNPs
        if set(list(haplotype)) == {"?"}:
            continue
        # remove single SNP
        if haplotype.count("?") == len(haplotype) - 1:
            continue

        compatible_count = 0
        compatible_seed = []
        hap_is_seed = haplotype in seeds

        for s in seeds:
            if is_compatible(haplotype, s):
                compatible_count += 1
                compatible_seed.append(s)
        
        if compatible_count == 0:
            cerr(f"Warning: haplotype {haplotype} could not be assigned to any seed, skipping.")
        else:
            new_count = round(count/compatible_count, 5)
            hap_info = {"haplotype": haplotype, "count": new_count}
            for s in compatible_seed:
                if s == haplotype:
                    # this hap is the seed of this group
                    # seed will be first and have the highest priority
                    groups[s]["h"] = [hap_info] + groups[s]["h"]
                elif hap_is_seed:
                    # this hap is used as seed of another group, least priority
                    groups[s]["t"].append(hap_info)
                else:
                    groups[s]["h"].append(hap_info)

    compatible_groups = []
    for s in seeds:
        final = groups[s]["h"] + groups[s]["t"]
        haplotypes = [hap["haplotype"] for hap in final]
        counts = [hap["count"] for hap in final]
        compatible_groups.append({"seed": s, "haplotypes": haplotypes, "counts": counts})
    
    return compatible_groups

def assemble_haplotype(df, max_consider = 15, prioritise = "count"):
    if prioritise not in ["count", "completeness"]:
        cerr(f"Invalid prioritise option: {prioritise}, must be one of ['count', 'n_missing']")
        cexit(1)
    
    subdf = df.copy(deep = True)
    subdf["n_missing"] = subdf["haplotype"].apply(lambda x: x.count("?"))

    if prioritise == "count":
        subdf.sort_values(by=["count", "n_missing"], ascending=[False, True], inplace=True)
    elif prioritise == "completeness":
        subdf.sort_values(by=["n_missing", "count"], ascending=[True, False], inplace=True)
    
    print("Scanning for seeds...")
    seeds = scan_for_seed(subdf["haplotype"].tolist())
    print(f"Found {len(seeds)} seeds")

    compatible_groups = assign_to_groups(subdf, seeds)
    
    assembled_haplotypes = []
    for group in compatible_groups:
        best_assembled = group["haplotypes"][0]
        support = group["counts"][0]
        haplotype_origins = [best_assembled]
        supporting_depth = [support if base != "?" else 0 for base in best_assembled]
        for haplotype, count in zip(group["haplotypes"][1:], group["counts"][1:]):
            new_assembled = []
            assembled_depth = [0 for _ in range(len(best_assembled))]
            for idx, (assembled_base, hap_base) in enumerate(zip(best_assembled, haplotype)):
                match assembled_base, hap_base:
                    case "?", _:
                        new_assembled.append(hap_base)
                        assembled_depth[idx] += count # add depth because 
                    case _, "?":
                        new_assembled.append(assembled_base)
                    case base1, base2 if base1 == base2:
                        new_assembled.append(base1)
                        assembled_depth[idx] += count
                    case base1, base2 if base1 != base2:
                        new_assembled = []
                        break
            if len(new_assembled) == 0:
                continue
            else:
                best_assembled = "".join(new_assembled)
                supporting_depth = [a + b for a, b in zip(supporting_depth, assembled_depth)]
                # support += count
                haplotype_origins.append(haplotype)
        assert len(best_assembled) == len(supporting_depth)
        assembled_haplotypes.append({"haplotype": best_assembled, "depths": supporting_depth, "assembled_from": haplotype_origins})
    return assembled_haplotypes

def main(args):
    import pandas as pd
    df = pd.read_csv(args.infile, sep="\t")
    
    if df["sample"].nunique() > 1:
        cexit("Input file contains multiple samples, please provide a single sample haplotype file.")
    
    if args.log is not None:
        log_fh = open(args.log, "w+")
        os.sys.stdout = log_fh
    
    haplotype_results = {}
    for marker in df["marker"].unique():
        print(f"Marker: {marker}")
        subdf = df[df["marker"] == marker]
        haplotype_results[marker] = assemble_haplotype(subdf, max_consider=args.max_consider, prioritise=args.prioritise)
    
    assembled_hap_df = pd.concat([
        pd.DataFrame(
            [(df["sample"].unique()[0] , marker, hap["haplotype"], hap["depths"], hap["assembled_from"]) for hap in haplotype_results[marker]],
            columns=["sample", "marker", "haplotype", "depths", "assembled_from"]
        ) for marker in haplotype_results
    ])

    assembled_hap_df.loc[:, "q25_depth"] = assembled_hap_df["depths"].apply(lambda r: pd.Series(r).quantile(0.25))
    assembled_hap_df[["sample", "marker", "haplotype", "q25_depth", "assembled_from", "depths"]] \
        .sort_values(["marker", "q25_depth", "haplotype"]).to_csv(args.outfile, sep="\t", index=False)
