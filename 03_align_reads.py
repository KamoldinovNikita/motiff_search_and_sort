#!/usr/bin/env python3
"""Выравнивание полных ридов."""

import argparse
import os
import shutil
import subprocess
from pathlib import Path


def run(args, log):
    source = args.input
    out = args.output
    _, index = prepare_reference(args.reference, out, log)
    files = input_files(source, "*.fastq.gz", "*.fastq")
    prefixes = [f.name.split(".")[0] for f in files]
    if len(set(prefixes)) != len(prefixes):
        raise ValueError("Input names have duplicate sample prefixes")
    for fastq, prefix in zip(files, prefixes):
        sam = out / f"{prefix}.sam"
        run_command(
            [
                "bowtie2",
                "-x",
                index,
                "-U",
                fastq,
                "-p",
                "4",
                "--very-sensitive-local",
                "-S",
                sam,
            ],
            log=log,
        )
        sort_bam(sam, out / f"{prefix}.sorted.bam", threads=4, log=log)
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


def prepare_reference(original, out, log):
    original = original.resolve()
    folder = out / "reference"
    folder.mkdir()
    reference = folder / original.name
    shutil.copy2(original, reference)
    index = original.with_suffix("")
    if not any(Path(str(index) + suffix).is_file() for suffix in (".1.bt2", ".1.bt2l")):
        index = reference.with_suffix("")
        run_command(["bowtie2-build", reference, index], log)
    return reference, index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
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
