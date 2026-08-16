# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

from __future__ import annotations

__copyright__ = "(C) 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

"""
generate_empty_gvcf.py

Generates an empty GVCF file suitable for GLnexus joint calling using pysam.
Accepts reference sequence metadata from either a FASTA index (.fai) file or
a SAM/BAM sequence dictionary (.dict) file.

Usage Examples:
    # Using --fai flag
    python3 generate_empty_gvcf.py --contig chrY -o output.g.vcf.gz --fai genome.fa.fai --sample-name SAMPLE1

    # Using --dict flag
    python3 generate_empty_gvcf.py --contig chrY -o output.g.vcf.gz --dict genome.dict --sample-name SAMPLE1

"""

import os
import sys
import argparse
import struct
import zlib
import typing
from ngs_pipeline import cerr, cexit, arg_parser, check_NGSENV_BASEDIR

if typing.TYPE_CHECKING:
    import pysam


def init_argparser():
    parser = arg_parser(
        desc="Generate empty GVCF files suitable for GLnexus joint calling.",
    )
    parser.add_argument(
        "--contig", type=str, help="Target contig name (e.g., chr1, chrY, 21)"
    )
    parser.add_argument(
        "--fai", type=str, help="Path to reference FASTA index (.fai) file"
    )
    parser.add_argument(
        "--dict",
        dest="dict_file",
        type=str,
        default=None,
        help="Path to SAM/BAM sequence dictionary (.dict) file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output filename (ends in .g.vcf.gz or .vcf.gz for BGZF compression)",
    )
    parser.add_argument(
        "-s",
        "--sample-name",
        type=str,
        default=None,
        help="Sample name to insert in #CHROM header line (default: derived from output filename)",
    )
    parser.add_argument(
        "-m",
        "--mode",
        choices=["header-only", "non-ref-block"],
        default="header-only",
        help="GVCF mode: 'header-only' (0 records) or 'non-ref-block' (includes 1..contig_length symbolic block)",
    )
    parser.add_argument(
        "--symbolic-allele",
        type=str,
        default="<NON_REF>",
        help="Symbolic non-reference allele used when --mode=non-ref-block (<NON_REF> or <*>)",
    )
    parser.add_argument(
        "--target-contig-only",
        action="store_true",
        help="Include only the target contig in VCF header rather than all contigs from .fai",
    )
    parser.add_argument(
        "--no-index",
        action="store_true",
        help="Disable automatic Tabix (.tbi) index file generation for compressed outputs",
    )

    return parser


# Maximum uncompressed block size for BGZF is 64 KB (65,536 bytes)
MAX_BGZF_BLOCK_SIZE = 60000


def compress_bgzf_block(data: bytes) -> bytes:
    """
    Compress a block of bytes (<= 64KB) into BGZF (Block Compressed GZIP) format.
    """
    if len(data) > 65280:
        raise ValueError(
            "Uncompressed block exceeds maximum BGZF capacity (65280 bytes)."
        )

    compressor = zlib.compressobj(level=6, method=zlib.DEFLATED, wbits=-15)
    compressed_data = compressor.compress(data) + compressor.flush()

    block_size = 18 + len(compressed_data) + 8
    bsize = block_size - 1

    # BGZF header (18 bytes)
    header = struct.pack(
        "<BBBBIBBHBBHH",
        31,
        139,
        8,
        4,  # ID1, ID2, CM (deflate), FLG (FEXTRA)
        0,
        0,
        255,  # MTIME, XFL, OS
        6,  # XLEN (subfield len)
        66,
        67,
        2,  # Subfield SI1='B', SI2='C', SLEN=2
        bsize,  # BSIZE (total block size - 1)
    )

    crc = zlib.crc32(data) & 0xFFFFFFFF
    isize = len(data) & 0xFFFFFFFF
    footer = struct.pack("<II", crc, isize)

    return header + compressed_data + footer


def get_bgzf_eof() -> bytes:
    """Return standard BGZF EOF marker block (28 bytes)."""
    return compress_bgzf_block(b"")


def create_tabix_index_for_empty_vcf(contigs: list[str]) -> bytes:
    """
    Generate a pure-Python Tabix (.tbi) index for an empty/header-only BGZF VCF file.
    """
    magic = b"TBI\1"
    n_ref = len(contigs)
    fmt = 2  # VCF format preset
    col_seq = 1  # #CHROM is column 1
    col_beg = 2  # POS is column 2
    col_end = 0  # 0 for single POS column
    meta = ord("#")  # Meta character
    skip = 0  # Header lines skip count

    hdr = struct.pack("<iIIIIii", n_ref, fmt, col_seq, col_beg, col_end, meta, skip)

    names_bytes = b"".join(c.encode("ascii") + b"\0" for c in contigs)
    l_nm = len(names_bytes)
    hdr += struct.pack("<i", l_nm) + names_bytes

    for _ in contigs:
        hdr += struct.pack("<ii", 0, 0)  # 0 bins, 0 intervals

    hdr += struct.pack("<Q", 0)  # n_no_coor = 0

    return compress_bgzf_block(hdr) + get_bgzf_eof()


def parse_fai_1(fai_path: str) -> dict[str, int]:
    """
    Parse a FASTA index (.fai) file.
    Returns an ordered dictionary mapping contig names to sequence lengths.
    """
    if not os.path.isfile(fai_path):
        raise FileNotFoundError(f"FASTA index (.fai) file not found at: '{fai_path}'")

    contigs = {}
    with open(fai_path, "r", encoding="utf-8") as f:
        for line_idx, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    f"Malformed .fai line {line_idx} in '{fai_path}': expected at least 2 tab-separated columns."
                )
            contig_name = parts[0]
            try:
                contig_len = int(parts[1])
            except ValueError:
                raise ValueError(
                    f"Invalid contig length on line {line_idx} of '{fai_path}': {parts[1]}"
                )
            contigs[contig_name] = contig_len

    if not contigs:
        raise ValueError(f".fai file '{fai_path}' contains no contigs.")

    return contigs


def parse_fai(fai_path: str) -> dict[str, int]:
    """
    Parse a FASTA index (.fai) file.
    Format: <contig_name>\t<length>\t<offset>\t<line_bases>\t<line_bytes>
    """
    if not os.path.isfile(fai_path):
        raise FileNotFoundError(f"FASTA index (.fai) file not found: '{fai_path}'")

    contigs = {}
    with open(fai_path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                raise ValueError(
                    f"Malformed .fai line {line_num} in '{fai_path}': expected at least 2 tab-delimited columns."
                )
            contig_name = cols[0]
            try:
                contig_len = int(cols[1])
            except ValueError:
                raise ValueError(
                    f"Invalid contig length on line {line_num} in '{fai_path}': {cols[1]}"
                )
            contigs[contig_name] = contig_len

    if not contigs:
        raise ValueError(f".fai file '{fai_path}' contains no contig entries.")
    return contigs


def parse_dict(dict_path: str) -> dict[str, int]:
    """
    Parse a SAM/BAM sequence dictionary (.dict) file (e.g. from Picard CreateSequenceDictionary or samtools dict).
    Format: @SQ\tSN:<contig_name>\tLN:<length>...
    """
    if not os.path.isfile(dict_path):
        raise FileNotFoundError(
            f"Sequence dictionary (.dict) file not found: '{dict_path}'"
        )

    contigs = {}
    with open(dict_path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or not line.startswith("@SQ"):
                continue
            fields = line.split("\t")
            sn = None
            ln = None
            for field in fields[1:]:
                if field.startswith("SN:"):
                    sn = field[3:]
                elif field.startswith("LN:"):
                    try:
                        ln = int(field[3:])
                    except ValueError:
                        raise ValueError(
                            f"Invalid length tag on line {line_num} in '{dict_path}': {field}"
                        )
            if sn is not None and ln is not None:
                contigs[sn] = ln
            else:
                raise ValueError(
                    f"Malformed @SQ header on line {line_num} in '{dict_path}': missing SN: or LN: tag."
                )

    if not contigs:
        raise ValueError(
            f"No @SQ sequence entries found in dictionary file: '{dict_path}'"
        )
    return contigs


def load_contigs(
    fai_path: str | None = None,
    dict_path: str | None = None,
    ref_path: str | None = None,
) -> dict[str, int]:
    """
    Load contigs and lengths from --fai, --dict, or auto-detected ref_path.
    """
    if fai_path:
        return parse_fai(fai_path)
    elif dict_path:
        return parse_dict(dict_path)
    elif ref_path:
        if ref_path.endswith(".dict") or ref_path.endswith(".dict.txt"):
            return parse_dict(ref_path)
        elif ref_path.endswith(".fai"):
            return parse_fai(ref_path)
        else:
            with open(ref_path, "r", encoding="utf-8") as f:
                first_line = f.readline()
            if first_line.startswith("@HD") or first_line.startswith("@SQ"):
                return parse_dict(ref_path)
            else:
                return parse_fai(ref_path)
    else:
        raise ValueError(
            "Must provide reference contig metadata via --fai, --dict, or positional argument."
        )


def build_gvcf_content(
    contig: str,
    fai_contigs: dict[str, int],
    sample_name: str,
    mode: str = "header-only",
    symbolic_allele: str = "<NON_REF>",
    target_contig_only: bool = False,
) -> str:
    """
    Construct VCF header and optional GVCF block record for GLnexus.
    """
    if contig and contig not in fai_contigs:
        available = ", ".join(list(fai_contigs.keys())[:10])
        if len(fai_contigs) > 10:
            available += f", ... ({len(fai_contigs) - 10} more)"
        raise ValueError(
            f"Contig '{contig}' was not found in the provided .fai file.\n"
            f"Available contigs in .fai: {available}"
        )

    lines = []
    lines.append("##fileformat=VCFv4.2")
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">')
    lines.append(
        '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele, for GVCF files">'
    )
    lines.append(
        '##ALT=<ID=<*>,Description="Represents any possible alternative allele, for GVCF files">'
    )

    # Standard INFO headers expected by GLnexus joint calling presets
    lines.append(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
    )
    lines.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')

    # Standard FORMAT headers expected by GLnexus joint calling presets
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append(
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
    )
    lines.append(
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">'
    )
    lines.append(
        '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP of positions merged in the block">'
    )
    lines.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">')
    lines.append(
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">'
    )
    lines.append(
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">'
    )
    lines.append(
        '##FORMAT=<ID=MED_DP,Number=1,Type=Integer,Description="Median DP of positions merged in the block">'
    )

    # Contig header definitions
    if contig:
        lines.append(f"##contig=<ID={contig},length={fai_contigs[contig]}>")
    else:
        for c_id, c_len in fai_contigs.items():
            lines.append(f"##contig=<ID={c_id},length={c_len}>")

    # Column header line
    lines.append(
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )

    # Optional non-ref symbolic block record
    if mode == "non-ref-block":
        contig_len = fai_contigs[contig]
        record = f"{contig}\t1\t.\tN\t{symbolic_allele}\t.\tPASS\tEND={contig_len}\tGT:DP:GQ:MIN_DP\t./.:0:0:0"
        lines.append(record)

    return "\n".join(lines) + "\n"


def build_pysam_header(
    fai_contigs: dict[str, int],
    sample_name: str,
    target_contig: str,
    target_contig_only: bool = False,
) -> pysam.VariantHeader:
    """
    Build a pysam.VariantHeader configured for GLnexus compatibility.
    """

    import pysam

    header = pysam.VariantHeader()
    header.add_sample(sample_name)

    # Standard filter and ALT definitions
    header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    header.add_line(
        '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele, for GVCF files">'
    )
    header.add_line(
        '##ALT=<ID=<*>,Description="Represents any possible alternative allele, for GVCF files">'
    )

    # Standard INFO headers expected by GLnexus
    header.add_line(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
    )
    header.add_line('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')

    # Standard FORMAT headers expected by GLnexus presets
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line(
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
    )
    header.add_line(
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">'
    )
    header.add_line(
        '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP of positions merged in the block">'
    )
    header.add_line(
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'
    )
    header.add_line(
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">'
    )
    header.add_line(
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">'
    )
    header.add_line(
        '##FORMAT=<ID=MED_DP,Number=1,Type=Integer,Description="Median DP of positions merged in the block">'
    )

    # Contig definitions
    if target_contig:
        header.add_line(
            f"##contig=<ID={target_contig},length={fai_contigs[target_contig]}>"
        )
    else:
        for c_id, c_len in fai_contigs.items():
            header.add_line(f"##contig=<ID={c_id},length={c_len}>")

    return header


def do_generate_empty_gvcf_1(
    contig: str,
    fai_path: str,
    output_filename: str,
    sample_name: str | None = None,
    mode: str = "header-only",
    symbolic_allele: str = "<NON_REF>",
    create_index: bool = True,
) -> None:
    """
    Main function to generate an empty GVCF file and optional Tabix index.
    """
    fai_contigs = parse_fai(fai_path)

    if not sample_name:
        base = os.path.basename(output_filename)
        for ext in [".g.vcf.gz", ".vcf.gz", ".g.vcf", ".vcf", ".gz", ".bgz"]:
            if base.endswith(ext):
                base = base[: -len(ext)]
                break
        sample_name = base if base else "EMPTY_SAMPLE"

    vcf_text = build_gvcf_content(
        contig=contig,
        fai_contigs=fai_contigs,
        sample_name=sample_name,
        mode=mode,
        symbolic_allele=symbolic_allele,
    )

    out_dir = os.path.dirname(os.path.abspath(output_filename))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    is_compressed = output_filename.endswith(".gz") or output_filename.endswith(".bgz")

    if is_compressed:
        raw_bytes = vcf_text.encode("utf-8")

        with open(output_filename, "wb") as f:
            offset = 0
            while offset < len(raw_bytes):
                chunk = raw_bytes[offset : offset + MAX_BGZF_BLOCK_SIZE]
                f.write(compress_bgzf_block(chunk))
                offset += len(chunk)
            f.write(get_bgzf_eof())

        print(f"Successfully generated BGZF compressed GVCF: {output_filename}")

        if create_index:
            tbi_filename = output_filename + ".tbi"
            contig_list = [contig] if contig else list(fai_contigs.keys())
            tbi_bytes = create_tabix_index_for_empty_vcf(contig_list)
            with open(tbi_filename, "wb") as f:
                f.write(tbi_bytes)
            print(f"Successfully generated Tabix index: {tbi_filename}")
    else:
        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(vcf_text)
        print(f"Successfully generated uncompressed GVCF: {output_filename}")


def do_generate_empty_gvcf(
    contig: str,
    output_filename: str,
    fai_path: str | None = None,
    dict_path: str | None = None,
    sample_name: str | None = None,
    mode: str = "header-only",
    symbolic_allele: str = "<NON_REF>",
    create_index: bool = True,
) -> None:
    """
    Generate an empty GVCF using pysam with reference metadata from .fai or .dict.
    """

    import pysam

    contigs_dict = load_contigs(fai_path=fai_path, dict_path=dict_path)

    if contig not in contigs_dict:
        available = ", ".join(list(contigs_dict.keys())[:10])
        if len(contigs_dict) > 10:
            available += f", ... ({len(contigs_dict) - 10} more)"
        raise ValueError(
            f"Contig '{contig}' was not found in reference dictionary/fai.\n"
            f"Available contigs: {available}"
        )

    if not sample_name:
        base = os.path.basename(output_filename)
        for ext in [".g.vcf.gz", ".vcf.gz", ".g.vcf", ".vcf", ".gz", ".bgz", ".bcf"]:
            if base.endswith(ext):
                base = base[: -len(ext)]
                break
        sample_name = base if base else "EMPTY_SAMPLE"

    header = build_pysam_header(
        fai_contigs=contigs_dict,
        sample_name=sample_name,
        target_contig=contig,
    )

    out_dir = os.path.dirname(os.path.abspath(output_filename))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    is_compressed = output_filename.endswith(".gz") or output_filename.endswith(".bgz")
    pysam_mode = "wz" if is_compressed else "w"

    with pysam.VariantFile(output_filename, mode=pysam_mode, header=header) as vcf_out:
        if mode == "non-ref-block":
            contig_len = contigs_dict[contig]
            # Create a symbolic non-ref block covering 1..contig_len
            record = vcf_out.new_record(
                contig=contig,
                start=0,  # 0-based coordinate in pysam (corresponds to POS 1 in VCF)
                stop=contig_len,
                alleles=("N", symbolic_allele),
                id=".",
                qual=None,
                filter="PASS",
                info={"END": contig_len},
            )
            # Genotype no-call with zero depth
            record.samples[sample_name]["GT"] = (None, None)
            record.samples[sample_name]["DP"] = 0
            record.samples[sample_name]["GQ"] = 0
            record.samples[sample_name]["MIN_DP"] = 0
            vcf_out.write(record)

    print(f"Successfully generated GVCF via pysam: {output_filename}")

    # Generate Tabix index for compressed files
    if is_compressed and create_index:
        try:
            pysam.tabix_index(output_filename, preset="vcf", force=True)
            print(
                f"Successfully generated Tabix index via pysam: {output_filename}.tbi"
            )
        except Exception:
            # Fallback for 0-record empty files if pysam tabix raises empty file notice
            tbi_filename = output_filename + ".tbi"
            contig_list = [contig] if target_contig_only else list(contigs_dict.keys())
            tbi_bytes = create_empty_tabix_index(contig_list)
            with open(tbi_filename, "wb") as f:
                f.write(tbi_bytes)
            print(f"Successfully generated Tabix index: {tbi_filename}")


def generate_empty_gvcf(args):

    print(args)

    try:
        do_generate_empty_gvcf(
            contig=args.contig,
            fai_path=args.fai,
            dict_path=args.dict_file,
            output_filename=args.output,
            sample_name=args.sample_name,
            mode=args.mode,
            symbolic_allele=args.symbolic_allele,
            create_index=not args.no_index,
        )
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


def create_spanning_empty_gvcf(contig_name, fai_path, output_vcf_path, sample):

    import pysam

    if not output_vcf_path.endswith(".gz"):
        output_vcf_path += ".gz"

    contig_length = None
    with open(fai_path, "r") as fai:
        for line in fai:
            parts = line.strip().split("\t")
            if parts[0] == contig_name:
                contig_length = int(parts[1])
                break

    if contig_length is None:
        print(f"Error: Contig '{contig_name}' not found.", file=sys.stderr)
        sys.exit(1)

    header = pysam.VariantHeader()
    header.version = "VCFv4.2"
    header.contigs.add(contig_name, length=contig_length)
    header.add_meta(
        "ALT",
        {"ID": "NON_REF", "Description": "Represents any possible alternative allele"},
    )  # type: ignore
    header.info.add("END", 1, "Integer", "Stop position of the interval")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("AD", "R", "Integer", "Allelic depths")
    header.formats.add("DP", 1, "Integer", "Read depth")
    header.formats.add("GQ", 1, "Integer", "Genotype Quality")
    header.formats.add("PL", "G", "Integer", "Phred-scaled likelihoods")
    header.add_sample(sample)

    with pysam.VariantFile(output_vcf_path, "wz", header=header) as vcf_out:
        # Create a single hom-ref block covering the entire contig
        record = vcf_out.new_record()
        record.chrom = contig_name
        record.pos = 1  # 1-based start
        record.stop = contig_length
        record.ref = "N"  # or reference base if known, 'N' acts as a placeholder
        record.alleles = ("N", "<NON_REF>")
        record.info["END"] = contig_length

        # Format fields for the sample
        record.samples[sample]["GT"] = (0, 0)
        record.samples[sample]["DP"] = 0
        record.samples[sample]["GQ"] = 0
        record.samples[sample]["AD"] = (0, 0)
        record.samples[sample]["PL"] = (0, 0, 0)

        vcf_out.write(record)

    pysam.tabix_index(output_vcf_path, preset="vcf", force=True)
    print(
        f"Successfully created spanning GVCF and index for {contig_name} (1-{contig_length})"
    )


def main(args):
    generate_empty_gvcf(args)


# EOF
