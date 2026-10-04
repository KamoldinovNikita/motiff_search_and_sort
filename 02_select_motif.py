#!/usr/bin/env python3
"""Отбор ридов с заданным мотивом."""

import argparse
import os
import subprocess
from pathlib import Path


MOTIF = "GCCAAGGC"


def run(args, log):
    source = args.input
    out = args.output
    for fastq in input_files(source, "*.fastq.gz"):
        target = out / f"{fastq.name[:-9]}.{args.motif}.fastq.gz"
        run_command(
            ["seqkit", "grep", "-s", "-P", "-p", args.motif, fastq, "-o", target],
            log=log,
        )
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--motif", default=MOTIF)
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
