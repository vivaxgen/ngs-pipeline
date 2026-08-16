# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

from __future__ import annotations

__copyright__ = "(C) 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
import typing

from ngs_pipeline import cerr, arg_parser

# command-line calculate-cpes-stats:


def init_argparser():
    p = arg_parser("collect stats from various steps")
    p.add_argument("-o", "--outfile", default="-", help="output filename [stdout]")
    p.add_argument("--trimmed", action="append")
    p.add_argument("--mapped", action="append")
    p.add_argument("--final", action="append")
    p.add_argument("--depth", action="append")
    p.add_argument("--mindepth", default=5, type=int, help="min depth to be stats [5]")
    p.add_argument("--sample", help="Sample code")
    p.add_argument("statfiles", nargs="*", help="stat files to be collected")
    return p


def parse_read_trimming(infiles):
    """return (no_of_original_reads, no_of_filtered_reads)"""

    # this will read a json file with the following specs:
    # {'original_reads': int, 'filtered_reads': int}

    import json

    original_reads = 0
    filtered_reads = 0

    for infile in infiles:
        with open(infile, "r") as f_in:
            d = json.load(f_in)
            original_reads += d["original_reads"]
            filtered_reads += d["filtered_reads"]

    return (original_reads, filtered_reads)


def parse_number(line, func: typing.Callable = int):
    # read number after colon
    return func(line.split(":")[-1].strip().split()[0].replace(",", ""))


class MappingStats(typing.NamedTuple):
    mapped: int
    proper: int
    insert_size: float
    insert_size_sd: float
    avg_qual: float
    inward: int
    outward: int
    other: int
    trans: int


def parse_mapping_stats(infiles):
    # read from output of samtools stats

    mapped = 0
    proper = 0
    insert_size = 0.0
    insert_size_sd = 0.0
    avg_qual = 0.0
    inward = 0
    outward = 0
    other = 0
    trans = 0

    N = len(infiles)

    for infile in infiles:
        cerr(f"[Reading stat file: {infile}]")
        with open(infile, "r") as f_in:
            for line in f_in:
                if not line.startswith("SN\t"):
                    continue
                if "reads mapped:" in line:
                    mapped += parse_number(line)
                    continue
                if "reads properly paired:" in line:
                    proper += parse_number(line)
                    continue
                if "insert size average:" in line:
                    insert_size += parse_number(line, float)
                    continue
                if "insert size standard deviation:" in line:
                    insert_size_sd += parse_number(line, float)
                    continue
                if "average quality:" in line:
                    avg_qual += parse_number(line, float)
                    continue
                if "inward oriented pairs:" in line:
                    inward += parse_number(line)
                    continue
                if "outward oriented pairs:" in line:
                    outward += parse_number(line)
                    continue
                if "pairs with other orientation:" in line:
                    other += parse_number(line)
                    continue
                if "pairs on different chromosomes:" in line:
                    trans += parse_number(line)
                    continue

    res = MappingStats(
        mapped=mapped,
        proper=proper,
        insert_size=insert_size / N,
        insert_size_sd=insert_size_sd / N,
        avg_qual=avg_qual / N,
        inward=inward * 2,
        outward=outward * 2,
        other=other * 2,
        trans=trans * 2,
    )

    return res


def calculate_depth(infile, mindepth=5):

    import gzip
    import statistics

    depths = []
    with gzip.open(infile, "r") as fin:
        cerr(f"[Reading depth file: {infile}]")
        for line in fin:
            d = int(line.strip().split()[-1])
            if d > mindepth:
                depths.append(d)

    total = sum(depths)
    L = len(depths)
    if total == 0:
        average = q1 = 0
    else:
        average = total / len(depths)
        q1 = statistics.quantiles(depths, n=4)[0]

    return (L, average, q1)


def calculate_depths(infiles, mindepth=5):

    import gzip
    import statistics

    depths = {}
    for infile in infiles:
        with gzip.open(infile, "r") as fin:
            cerr(f"[Reading depth file: {infile}]")
            for line in fin:
                tokens = line.strip().split()
                try:
                    chr_d = depths[tokens[0]]
                except KeyError:
                    chr_d = depths[tokens[0]] = {}

                d = int(tokens[-1])
                try:
                    chr_d[tokens[1]] += d
                except KeyError:
                    chr_d[tokens[1]] = d

    depth_values = []

    for chr_d in depths.values():
        for d in chr_d.values():
            if d >= mindepth:
                depth_values.append(d)

    total = sum(depth_values)
    L = len(depth_values)
    if total == 0:
        average = q1 = 0
    else:
        average = total / L
        q1 = statistics.quantiles(depth_values, n=4)[0]

    return (L, average, q1)


def calculate_cpes_stats(args):

    initial_reads, trimmed_reads = parse_read_trimming(args.trimmed)

    # if initial_reads is zero, set it to 1 to avoid division by zero
    initial_reads = max(1, initial_reads)
    trimmed_r = trimmed_reads / initial_reads

    # avoid division by zero for continuing calculations
    trimmed_reads = max(1, trimmed_reads)

    # get initial mapping stats
    mapping_stat = parse_mapping_stats(args.mapped)

    # for each input statfiles, make it as filename::stage-{idx}

    statfile_by_len = {}

    for statfile in args.statfiles:
        fn, stage = statfile.split("::")
        stage, idx = stage.rsplit("-", 1)
        stage = stage.removesuffix("_")

        fn_len = len(fn)
        if fn_len not in statfile_by_len:
            statfile_by_len[fn_len] = (stage, [fn])
        else:
            if statfile_by_len[fn_len][0] != stage:
                raise ValueError(
                    f"Stage mismatch for files with same length: {statfile_by_len[fn_len][0]} vs {stage}"
                )
            statfile_by_len[fn_len][1].append(fn)

    results_by_len = {}
    for k, v in statfile_by_len.items():
        stage, fns = v
        cerr(f"[Stage: {stage}, files: {len(fns)}]")
        results_by_len[k] = (stage, parse_mapping_stats(fns))

    sorted_results = sorted(results_by_len.items(), key=lambda x: x[0])

    headers = [
        "SAMPLE",
        "RAW",
        "TRIMMED",
        "TRIMMED_R",
        "MAPPED",
        "MAPPED_R",
    ]

    values = [
        args.sample,
        f"{initial_reads}",
        f"{trimmed_reads}",
        f"{trimmed_r:5.3f}",
        f"{mapping_stat.mapped}",
        f"{mapping_stat.mapped / trimmed_reads:5.3f}",
    ]

    prev_stats = mapping_stat
    prev_stage = "mapped"
    for k, v in sorted_results:
        stage, stats = v

        if stage.startswith("filter"):
            # get the orientation stats of the previous state
            headers.extend(
                [
                    f"{prev_stage}_INWARD",
                    f"{prev_stage}_INWARD_R",
                    f"{prev_stage}_OUTWARD",
                    f"{prev_stage}_OUTWARD_R",
                    f"{prev_stage}_OTHER",
                    f"{prev_stage}_OTHER_R",
                    f"{prev_stage}_TRANS",
                    f"{prev_stage}_TRANS_R",
                ]
            )
            values.extend(
                [
                    f"{prev_stats.inward}",
                    (
                        f"{prev_stats.inward / prev_stats.mapped:5.3f}"
                        if prev_stats.mapped > 0
                        else "0.000"
                    ),
                    f"{prev_stats.outward}",
                    (
                        f"{prev_stats.outward / prev_stats.mapped:5.3f}"
                        if prev_stats.mapped > 0
                        else "0.000"
                    ),
                    f"{prev_stats.other}",
                    (
                        f"{prev_stats.other / prev_stats.mapped:5.3f}"
                        if prev_stats.mapped > 0
                        else "0.000"
                    ),
                    f"{prev_stats.trans}",
                    (
                        f"{prev_stats.trans / prev_stats.mapped:5.3f}"
                        if prev_stats.mapped > 0
                        else "0.000"
                    ),
                ]
            )

        headers.extend(
            [
                f"{stage}",
                f"{stage}_R",
            ]
        )
        values.extend(
            [
                f"{stats.mapped}",
                (
                    f"{stats.mapped / prev_stats.mapped:5.3f}"
                    if prev_stats.mapped > 0
                    else "0.000"
                ),
            ]
        )

        prev_stats = stats
        prev_stage = stage

    # get the final stats

    basepairs, avg_depth, q1_depth = calculate_depths(args.depth, args.mindepth)

    headers.extend(
        [
            "FINAL_R",
            "INSERT_SIZE",
            "INSERT_SIZE_SD",
            "AVG_QUAL",
            "BASEPAIRS",
            "AVG_DEPTH",
            "Q1_DEPTH",
        ]
    )
    values.extend(
        [
            f"{prev_stats.mapped / trimmed_reads:5.3f}",
            f"{prev_stats.insert_size:5.3f}",
            f"{prev_stats.insert_size_sd:5.3f}",
            f"{prev_stats.avg_qual:5.3f}",
            f"{basepairs}",
            f"{avg_depth:5.2f}",
            f"{q1_depth:5.2f}",
        ]
    )

    # save to file
    outfile = sys.stdout if args.outfile == "-" else open(args.outfile, "w")
    try:
        print("\t".join(headers), file=outfile)
        print("\t".join(values), file=outfile)
    finally:
        if outfile is not sys.stdout:
            outfile.close()


def main(args):
    calculate_cpes_stats(args)


# EOF
