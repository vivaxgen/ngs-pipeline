# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

from __future__ import annotations

__copyright__ = "(C) 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"


import argparse
import asyncio
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

from ngs_pipeline import check_VVG_BASEDIR, cerr, cexit

# Type aliases (Python >=3.12)

type SampleHeaders = dict[int, str]
type HeaderDB = dict[str, SampleHeaders]


# Configuration
@dataclass(frozen=True, slots=True)
class Config:
    reads_subdir: str = "reads"
    input_prefix: str = "raw"
    output_prefix: str = "model"
    max_concurrency: int = 64
    sample_numbers: int = -1
    file_numbers: int = -1

    raw_pattern: re.Pattern[str] = field(init=False)

    def __post_init__(self):
        object.__setattr__(
            self,
            "raw_pattern",
            re.compile(rf"{re.escape(self.input_prefix)}-(\d+)\.fastq\.gz$"),
        )

    @property
    def input_glob(self) -> str:
        return f"{self.input_prefix}-*.fastq.gz"


# Data model
@dataclass(slots=True)
class Sample:
    name: str
    path: Path
    reads_dir: Path
    raw_files: list[Path]


# Discovery
def discover_samples(
    input_dirs: list[Path],
    config: Config,
) -> list[Sample]:
    samples: list[Sample] = []

    for root in input_dirs:
        if not root.is_dir():
            raise NotADirectoryError(root)

        for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            reads_dir = sample_dir / config.reads_subdir

            if not reads_dir.is_dir():
                continue

            raw_files = sorted(reads_dir.glob(config.input_glob))

            if not raw_files:
                continue

            samples.append(
                Sample(
                    name=f"{root.name}/{sample_dir.name}",
                    path=sample_dir,
                    reads_dir=reads_dir,
                    raw_files=raw_files,
                )
            )

    return samples


# Processing

ONT_MODEL_re = re.compile(r"RG:Z:[a-f0-9-]+_((?:dna|rna)_[a-z0-9._]+@v[0-9.]+)")


def process_header(line: str) -> str:
    """return the ONT model from header line, or empty string if not found"""
    match = ONT_MODEL_re.search(line)
    if match:
        return (
            match.group(1)
            .removeprefix("dna_")
            .removeprefix("rna_")
            .replace(".", "")
            .replace("@", "_")
        )
    return ""


def process_dictionary(headers: HeaderDB) -> int:
    """get the set of the ONT models and check if they are already present."""

    # collect the models from the headers
    N = 0
    models: set[str] = set()
    for sample_headers in headers.values():
        for header in sample_headers.values():
            if header:
                models.add(header)
                N += 1

    cerr(f"Found ONT models: {sorted(models)}")

    vvg_basedir = check_VVG_BASEDIR()

    missing_models = []
    base_dir = Path(vvg_basedir) / "opt" / "clair3-models"
    for model in models:
        model_dir = base_dir / model
        cerr(f"Checking for model '{model}' in {model_dir}...")
        if not model_dir.is_dir():
            missing_models.append(model)

    if any(missing_models):
        cerr(f"Missing ONT models: {sorted(missing_models)}")
        cexit(
            "Please fetch the missing models using the 'fetch_clair3_models' command."
        )

    return N


# Blocking helpers


def read_first_line(path: Path) -> str:
    with gzip.open(path, "rt") as fp:
        return fp.readline()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w") as fp:
        fp.write(text)
        fp.write("\n")


# Async helpers


async def read_one_file(
    path: Path,
    config: Config,
    sem: asyncio.Semaphore,
) -> tuple[int, str]:
    m = config.raw_pattern.match(path.name)

    if m is None:
        raise ValueError(f"Unexpected filename: {path}")

    idx = int(m.group(1))

    async with sem:
        line = await asyncio.to_thread(
            read_first_line,
            path,
        )

    return idx, process_header(line)


async def write_one_file(
    sample: Sample,
    idx: int,
    header: str,
    config: Config,
    sem: asyncio.Semaphore,
) -> None:
    out = sample.reads_dir / f"{config.output_prefix}-{idx}.txt"

    async with sem:
        await asyncio.to_thread(
            write_text,
            out,
            header,
        )


# Sample tasks


async def read_sample(
    sample: Sample,
    config: Config,
    sem: asyncio.Semaphore,
) -> tuple[str, SampleHeaders]:
    results = await asyncio.gather(
        *(
            read_one_file(
                raw,
                config,
                sem,
            )
            for raw in sample.raw_files
        )
    )

    return sample.name, dict(results)


async def write_sample(
    sample: Sample,
    headers: SampleHeaders,
    config: Config,
    sem: asyncio.Semaphore,
) -> None:
    await asyncio.gather(
        *(
            write_one_file(
                sample,
                idx,
                header,
                config,
                sem,
            )
            for idx, header in headers.items()
        )
    )


# Pipeline


async def run(
    input_dirs: list[Path],
    config: Config,
) -> None:
    samples = discover_samples(
        input_dirs,
        config,
    )

    if not samples:
        cexit("ERR:No samples found.")

    if config.sample_numbers > 0 and len(samples) != config.sample_numbers:
        cexit(
            f"ERR: Expected {config.sample_numbers} samples, "
            f"but found {len(samples)} samples."
        )

    sem = asyncio.Semaphore(config.max_concurrency)

    cerr(f"Found {len(samples)} samples")

    # Stage 1
    results = await asyncio.gather(
        *(
            read_sample(
                sample,
                config,
                sem,
            )
            for sample in samples
        )
    )

    headers: HeaderDB = dict(results)

    # Stage 2
    no_of_files = process_dictionary(headers)
    if config.file_numbers > 0 and no_of_files != config.file_numbers:
        cexit(
            f"ERR: Expected {config.file_numbers} files, "
            f"but processing {no_of_files} files."
        )

    # Stage 3
    await asyncio.gather(
        *(
            write_sample(
                sample,
                headers[sample.name],
                config,
                sem,
            )
            for sample in samples
        )
    )

    cerr("Done")


# CLI


def init_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--reads-subdir",
        default="reads",
    )

    parser.add_argument(
        "--input-prefix",
        default="raw",
    )

    parser.add_argument(
        "--output-prefix",
        default="model",
    )

    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=64,
    )

    parser.add_argument(
        "-s",
        "--sample-numbers",
        type=int,
        default=-1,
        help="number of samples for checking purpose (default: -1, meaning no checking)",
    )

    parser.add_argument(
        "-n",
        "--file-numbers",
        type=int,
        default=-1,
        help="number of files for checking purpose (default: -1, meaning no checking)",
    )

    parser.add_argument(
        "input_dirs",
        nargs="+",
        type=Path,
    )

    return parser


def check_clair3_models(args) -> None:

    config = Config(
        reads_subdir=args.reads_subdir,
        input_prefix=args.input_prefix,
        output_prefix=args.output_prefix,
        max_concurrency=args.jobs,
        sample_numbers=args.sample_numbers,
    )

    # run asyncio event loop
    asyncio.run(
        run(
            input_dirs=args.input_dirs,
            config=config,
        )
    )


def main(args) -> None:
    check_clair3_models(args)


# EOF
