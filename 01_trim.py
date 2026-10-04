#!/usr/bin/env python3
"""Тримминг FASTQ."""

import argparse
import json
import os
import subprocess
from pathlib import Path


def run(args, log):
    source = args.input
    out = args.output
    files = input_files(source, "*.fastq", "*.fastq.gz")
    prefixes = [f.name.split(".")[0] for f in files]
    if len(set(prefixes)) != len(prefixes):
        raise ValueError("Input names have duplicate prefixes before the first dot")
    with open(out / "trimming_summary.tsv", "w") as summary:
        summary.write("sample\treads_before\treads_after\tpercent_retained\n")
        for fastq, prefix in zip(files, prefixes):
            report = out / f"{prefix}.fastp.json"
            run_command(
                [
                    "fastp",
                    "-i",
                    fastq,
                    "-o",
                    out / f"{prefix}.trimmed.fastq.gz",
                    "-j",
                    report,
                    "-h",
                    out / f"{prefix}.fastp.html",
                    "--thread",
                    "4",
                ],
                log=log,
            )
            stats = json.loads(report.read_text())["summary"]
            before = stats["before_filtering"]["total_reads"]
            after = stats["after_filtering"]["total_reads"]
            percent = f"{after / before * 100:.2f}" if before else "0.00"
            summary.write(f"{prefix}\t{before}\t{after}\t{percent}\n")
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
