#!/usr/bin/env python3
"""Исключение 45 интервалов."""

import argparse
import os
import subprocess
from pathlib import Path


REGION_STARTS = (
    11899,
    296538,
    380975,
    426425,
    455160,
    476510,
    547352,
    741498,
    766201,
    891650,
    946130,
    1174563,
    1288329,
    1292585,
    1336490,
    1362606,
    1364014,
    1564469,
    1824170,
    2385537,
    2414808,
    2478061,
    2534210,
    2626436,
    2631312,
    3059099,
    3068402,
    3235449,
    3442625,
    3717120,
    3732665,
    3930148,
    3994789,
    4012687,
    4059244,
    4091613,
    4181470,
    4233270,
    4273026,
    4305851,
    4397170,
    4457743,
    4556104,
    4590077,
    4626345,
)


def run(args, log):
    source = args.input
    out = args.output
    bed = out / "regions_to_remove.bed"
    with open(bed, "w") as handle:
        for start in REGION_STARTS:
            handle.write(f"Ecoli_MG1655\t{start}\t{start + 8}\n")
    for bam in input_files(source, "*.bam"):
        temporary = out / f"{bam.stem}.kept.bam"
        command = ["bedtools", "intersect", "-abam", bam, "-b", bed, "-v"]
        with open(temporary, "wb") as handle:
            subprocess.run(command, stdout=handle, stderr=log, check=True)
        sort_bam(temporary, out / f"{bam.stem}.no_regions.sorted.bam", log=log)
        temporary.unlink()
    return out


def input_files(source, *patterns):
    if source.is_file():
        files = [source] if any(source.match(pattern) for pattern in patterns) else []
    else:
        files = []
        for pattern in patterns:
            files.extend(sorted(source.glob(pattern)))
    if not files:
        raise FileNotFoundError(f"No input files matching {patterns} in {source}")
    return files


def run_command(command, log):
    subprocess.run(command, stdout=log, stderr=log, check=True)


def sort_bam(sam, bam, log, threads=None):
    command = ["samtools", "sort"]
    if threads:
        command += ["-@", str(threads)]
    run_command(command + ["-o", bam, sam], log)
    run_command(["samtools", "index", bam], log)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bin-dir", action="append", type=Path, default=[])
    args = parser.parse_args()
    os.environ["PATH"] = os.pathsep.join(
        [str(path.resolve()) for path in args.bin_dir] + [os.environ["PATH"]]
    )
    args.input = args.input.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True)
    with open(args.output / "commands.log", "a") as log:
        run(args, log)


if __name__ == "__main__":
    main()
