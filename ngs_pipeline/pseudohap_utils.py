# Utility functions for pseudohap tables

import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings

NON_SAMPLE_COLS = ["marker", "haplotypes", "total_count", "num_samples", "presence_ratio", "presence_ratio_relative", "n_haplotype", "total_depth", "presence_n_sample"]

def marker_missingness(wide_sample_df, sample_cols=None):
    if sample_cols is None:
        sample_cols = set(wide_sample_df.columns).difference(NON_SAMPLE_COLS)
    missing_counts = wide_sample_df[list(sample_cols)].map(lambda x: 1 if any([y > 0 for y in map(int, str(x).split(",")) ]) else 0 ).apply(lambda x: sum([1 if y == 0 else 0 for y in x]), axis =1 )
    return pd.DataFrame({
        "marker": wide_sample_df["marker"],
        "missing_count": missing_counts,
        "missing_ratio": np.round(missing_counts/len(sample_cols), 2)
    })

def wide_hap_remove_empty_haplotypes(wide_df):
    if "sample" in wide_df.columns:
        wide_df = wide_df.set_index("sample")
    hap_to_drop = (wide_df.map(lambda x: 1 if x > 0 else 0).sum(axis =0).to_frame("count").query("count == 0")).index
    result_df = wide_df.drop(columns=hap_to_drop)
    warnings.warn(f"Dropping {len(hap_to_drop)} haplotypes that are not present in any sample.")
    return(result_df)

def wide_sample_to_wide_haplotype(wide_df, progress=False, sample_cols=None):
    if sample_cols is None:
        sample_cols = [col for col in wide_df.columns if not col in NON_SAMPLE_COLS]
    sample_data = {sample: {} for sample in sample_cols}
    for i, row in tqdm(wide_df.iterrows(), total=len(wide_df), disable=progress):
        haplotypes = row["haplotypes"].split(",")
        marker = row["marker"]
        new_df_cols = [f"{marker}|{hap}" for hap in haplotypes]
        for sample in sample_cols:
            val = row[sample]
            if pd.isna(val):
                counts = [0] * len(haplotypes)
            else:
                counts = [int(c) for c in str(val).split(",")]

            sample_data[sample].update(dict(zip(new_df_cols, counts)))
    result_df = pd.DataFrame.from_dict(sample_data, orient="index")
    result_df.index.name = "sample"
    result_df = result_df.reset_index()
    hap_to_drop = result_df.drop(columns=["sample"]).apply(lambda x: (x > 0).sum(), axis=0).to_frame(name="n_sample").query("n_sample < 1").index.tolist()
    result_df = result_df.drop(columns=hap_to_drop)
    warnings.warn(f"Dropping {len(hap_to_drop)} haplotypes that are not present in any sample.")
    return(result_df)

def supplement_wide_sample_df(wide_df, n_sig=2, sample_cols=None):
    results = []
    if sample_cols is None:
        sample_rows = [col for col in wide_df.columns if not col in NON_SAMPLE_COLS ]
    else:
        sample_rows = sample_cols
    n_sample = len(sample_rows)
    for i, row in wide_df.iterrows():
        hap_sample_counts = row[sample_rows].str.split(",")
        presence_sample = np.array(hap_sample_counts.apply(lambda x: [1 if int(a) > 0 else 0 for a in x]).tolist()).sum(axis = 0)
        total_depth = np.array(hap_sample_counts.apply(lambda x: [int(a) for a in x]).tolist()).sum(axis = 0)
        new_row = row.copy()
        new_row["num_samples"] = ",".join(presence_sample.astype(str))
        new_row["presence_ratio"] = ",".join(np.round(presence_sample / n_sample, decimals=n_sig).astype(str))
        new_row["presence_ratio_relative"] = ",".join(np.round(presence_sample / np.sum(presence_sample), decimals=n_sig).astype(str))
        new_row["total_count"] = ",".join(total_depth.astype(str))
        results.append(new_row)
    return pd.concat(results, axis=1).T[["marker", "haplotypes", "num_samples", "presence_ratio", "presence_ratio_relative", "total_count"] + sample_rows]

def wide_sample_to_major_wide_sample(wide_df, sample_cols=None):
    if sample_cols is None:
        sample_cols = [col for col in wide_df.columns if not col in NON_SAMPLE_COLS]

    def _mask_non_major(row):
        result = {}
        for col, cell in row.items():
            if not col in sample_cols:
                result[col] = cell
                continue
            hap_counts = np.array([int(v) for v in cell.split(",")])
            minor_idx = np.argsort(hap_counts, kind="stable")[:-1]
            hap_counts[minor_idx] = 0
            result[col] = ",".join(map(str, hap_counts))
        return pd.Series(result)

    result_df = wide_df.apply(_mask_non_major, axis = 1)

    return(result_df)

def trim_haplotypes_wide_sample(wide_df, min_count = 0, min_ratio=0, sample_cols=None, use_relative_ratio=False):
    if not set(["num_samples", "presence_ratio"]).issubset(wide_df.columns):
        new_df = supplement_wide_sample_df(wide_df, sample_cols=sample_cols)
    else:
        new_df = wide_df.copy()

    if sample_cols is None:
        sample_cols = [col for col in wide_df.columns if not col in NON_SAMPLE_COLS]

    allowed_non_samples_cols = new_df.columns.difference(sample_cols)

    def _trim_haplotype(row):
        haplotype = row["haplotypes"].split(",")
        num_samples = map(int, str(row["num_samples"]).split(","))
        if use_relative_ratio:
            presence_ratio = map(float, str(row["presence_ratio_relative"]).split(","))
        else:
            presence_ratio = map(float, str(row["presence_ratio"]).split(","))
        sample_cols = [c for c in row.index if not c in allowed_non_samples_cols]
        passed = [
            (i, hap, presence, ratio)
            for i, (hap, presence, ratio) in enumerate(zip(haplotype, num_samples, presence_ratio))
            if presence >= min_count and ratio >= min_ratio
        ]
        if not passed:
            passed_i, passed_haps, passed_presences, passed_ratios = [], [], [], []
        else:
            passed_i, passed_haps, passed_presences, passed_ratios = zip(*passed)
        if len(passed_haps) == 0:
            return None
        row["haplotypes"] = ",".join(map(str, passed_haps))
        row["num_samples"] = ",".join(map(str, passed_presences))
        row["presence_ratio"] = ",".join(map(str, passed_ratios))
        row["presence_ratio_relative"] = ",".join(map(str, [float(ratio) / sum(passed_ratios) if sum(passed_ratios) > 0 else 0 for ratio in passed_ratios]))
        for s in sample_cols:
            sample_vals = str(row[s]).split(",")
            filtered_vals = [sample_vals[i] for i in passed_i]
            row[s] = ",".join(map(str, filtered_vals))
        return row

    new_df = new_df.apply(_trim_haplotype, axis = 1).dropna(subset=['haplotypes']).reset_index(drop=True)
    return new_df

def long_to_wide_ordered(long_df, progress=False):
    # haplotypes ordered from left to right based on how common they are, followed by count, <-> lexical break tie
    passed_haplotype = long_df.query("passed_filters == True")
    samples_order = long_df["sample"].unique().tolist()

    marker_stats = long_df.groupby(["marker", "haplotype",  "passed_filters"]).agg(
            num_samples = ("sample", "nunique"),
            total_count = ("count", "sum"),
            min_depth = ("min_depth", "min")
        ).reset_index()

    sample_stats = long_df.query("passed_filters == True").groupby("sample").agg(
            num_markers = ("marker", "nunique")
        ).reset_index()

    sorted_marker = marker_stats.sort_values(["marker", "passed_filters", "num_samples", "total_count", "min_depth", "haplotype"],
            ascending=[True, False, False, False, False, True])[["marker", "haplotype", "passed_filters", "total_count", "num_samples"]]

    rows = []
    for marker in tqdm(sorted_marker["marker"].unique(), disable=progress):
        marker_df = sorted_marker.query("marker == @marker and passed_filters == True")
        haplotypes = marker_df["haplotype"].to_list()
        if haplotypes == []:
            haplotypes = ["*"]
        total_counts = marker_df["total_count"].astype(str).to_list()
        if total_counts == []:
            total_counts = ["0"]
        num_samples = marker_df["num_samples"].astype(str).to_list()
        if num_samples == []:
            num_samples = ["0"]
        marker_info = pd.DataFrame({"marker": marker, "haplotypes": ",".join(haplotypes), \
            "total_count": ",".join(total_counts), \
            "num_samples": ",".join(num_samples)}, index= [marker])
        long_df_marker = long_df.query("passed_filters == True and marker == @marker")
        if long_df_marker.empty:
            long_df_marker = pd.DataFrame({
                "sample": samples_order,
                "haplotype": ["*"] * len(samples_order),
                "count": [0] * len(samples_order)
            })
        pivoted = long_df_marker.pivot(index="haplotype", columns = "sample", values = "count") \
            .fillna(0).astype(int)
        sample_missing_cols = list(set(samples_order) - set(pivoted.columns.tolist()))
        pivoted.loc[:, sample_missing_cols] = 0

        pivoted = pivoted.loc[haplotypes, samples_order]
        pivoted = pivoted.apply(lambda x: ",".join(map(str, x)), axis = 0).to_frame(name=marker).T
        pivoted.columns.name = None
        rows.append(pd.concat([marker_info, pivoted], axis=1))
    wide_tab = pd.concat(rows, axis=0)
    return wide_tab

def process_individual_sample_haplotype(df, marker_df, min_depth_per_base = 5, haplotype_pass_min_count = 10, haplotype_pass_min_ratio = 0.1, sample_id="sample", use_poisson_threshold=False, error_rate=0.1, confidence=0.99):
    # sample, marker, haplotype, q25_depth, assembled_from, depths
    missing_markers = set(marker_df["marker"].unique()) - set(df["marker"].unique())

    if df["sample"].nunique() > 1:
        raise ValueError(f"Expected data for a single sample, but found {df['sample'].nunique()} samples in the input DataFrame.")
    if df.shape[0] > 0 and df["sample"].unique()[0] != sample_id:
        raise ValueError(f"Expected sample ID '{sample_id}', but found '{df['sample'].unique()[0]}' in the input DataFrame.")

    if len(missing_markers) > 0:
        df_base = []
        for marker in missing_markers:
            base_df = pd.DataFrame({
                "sample": sample_id,
                "marker": marker,
                "haplotype": ["?" * marker_df.query("marker == @marker").shape[0]],
                "q25_depth": [0],
                "assembled_from": [[]],
                "depths": [str([0 for _ in range(marker_df.query("marker == @marker").shape[0])])]
            })
            df_base.append(base_df)
        df = pd.concat([df] + df_base, ignore_index=True)
    final_common_cols = ["sample", "marker", "haplotype", "count", "min_depth", "passed_filters", "failed_reason"]
    failed_encoding = {
        0: "",
        1: "Support",
        2: "Incomplete",
        4: "Ratio",
        3: "Support;Incomplete",
        5: "Support;Ratio",
        6: "Incomplete;Ratio",
        7: "Support;Incomplete;Ratio"
    }
    df["passed_filters"] = True
    df["failed_reason"] = 0

    assert "q25_depth" in df.columns, "Expected 'q25_depth' column in assembled haplotype data."
    df["min_depth"] = df.apply(lambda row: min(eval(row["depths"])) if pd.notna(row["depths"]) else 0, axis=1)
    df["count"] = df["q25_depth"].astype(int)

    filtered_df = df[final_common_cols].copy()
    # 1. filter incomplete haplotype
    filtered_df.loc[:, "incomplete_haplotype"] = filtered_df["haplotype"].apply(lambda x: 'N' in x or '?' in x)
    filtered_df.loc[filtered_df["incomplete_haplotype"].values, "passed_filters"] = False
    filtered_df.loc[filtered_df["incomplete_haplotype"].values, "failed_reason"] += 2

    # 2. filter haplotype with low support (min depth per base < min_depth_per_base)
    filtered_df.loc[filtered_df.query("min_depth < @min_depth_per_base").index, "passed_filters"] = False
    filtered_df.loc[filtered_df.query("min_depth < @min_depth_per_base").index, "failed_reason"] += 1

    def get_poisson_threshold(counts, min_absolute_count=haplotype_pass_min_count, error_rate=haplotype_pass_min_ratio, confidence=0.99):
        import scipy.stats as st
        if all(count == 0 for count in counts):
            return 0
        major_count = max(counts)
        # Calculate expected noise reads generated by the major haplotype
        expected_noise = major_count * error_rate
        # Calculate the statistical cutoff
        dynamic_cutoff = st.poisson.ppf(confidence, expected_noise)
        # Return whichever is higher: the absolute minimum or the dynamic cutoff
        return max(min_absolute_count, dynamic_cutoff)

    # per marker total count max (haplotype_pass_min_count, haplotype_pass_min_ratio * total count of haplotypes that passed support and completeness filters)
    min_per_marker_ratio = np.ceil((filtered_df.query("passed_filters == True").groupby("marker")["count"].max() * haplotype_pass_min_ratio).clip(lower = haplotype_pass_min_count)).astype(int)

    if use_poisson_threshold:
        min_per_marker_poisson = filtered_df.query("passed_filters == True").groupby("marker")["count"].apply(lambda x: get_poisson_threshold(x, haplotype_pass_min_count, error_rate, confidence)).astype(int)
        min_per_marker = pd.Series(np.max([min_per_marker_poisson, min_per_marker_ratio], axis=0), index=min_per_marker_poisson.index)
    else:
        min_per_marker = min_per_marker_ratio

    filtered_df = filtered_df.merge(min_per_marker.rename("min_count"), on="marker")

    # 3. filter haplotype that do not meet the minimum count requirement (either absolute or relative)
    filtered_df.loc[filtered_df.query("count < min_count").index, "passed_filters"] = False
    filtered_df.loc[filtered_df.query("count < min_count").index, "failed_reason"] += 4

    stats = {"sample": sample_id,
        "passed": filtered_df.query("passed_filters == True").groupby("marker").size().to_dict()}

    filtered_df["failed_reason"] = filtered_df['failed_reason'].map(failed_encoding)
    # filtered_df[final_common_cols].to_csv(output[0], sep="\t", index=False)
    stats["n_markers_with_passed_haplotypes"] = (filtered_df.query("passed_filters == True").groupby("marker").size() > 0).sum()
    stats["naive_coi"] = filtered_df.query("passed_filters == True").groupby("marker").size().max()
    return {"df": filtered_df[final_common_cols], "stat": stats}

def generate_n_haplotype_stats(before_df, after_df=None):
    before_stats = before_df[["marker", "haplotypes"]].copy()
    before_stats["n_haplotype"] = before_stats["haplotypes"].apply(lambda x: len(x.split(",")) if pd.notna(x) else 0)
    if after_df is None:
        return before_stats[["marker", "n_haplotype"]]
    after_stats = after_df[["marker", "haplotypes"]].copy()
    after_stats["n_haplotype"] = after_stats["haplotypes"].apply(lambda x: len(x.split(",")) if pd.notna(x) else 0)
    stats_df = pd.DataFrame({
        "marker": before_stats["marker"],
        "n_haplotype_before": before_stats["n_haplotype"],
        "n_haplotype_after": after_stats["n_haplotype"]
    })
    return stats_df