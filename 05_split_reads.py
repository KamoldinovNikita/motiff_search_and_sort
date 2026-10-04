#!/usr/bin/env python3
"""Разделение start, middle20, end."""

import argparse
import os
import subprocess
from pathlib import Path


MOTIF = "GCCAAGGC"

K = 20


def run(args, log):
    for bam in input_files(args.input, "*.bam"):
        paths = {
            tag: args.output / f"{bam.stem}.{tag}.sam"
            for tag in ("start", "middle20", "end")
        }
        with (
            open(paths["start"], "w") as start,
            open(paths["middle20"], "w") as middle,
            open(paths["end"], "w") as end,
        ):
            files = {"start": start, "middle20": middle, "end": end}
            for line in read_lines(["samtools", "view", "-h", bam], log):
                if line.startswith("@"):
                    for file in files.values():
                        file.write(line)
                    continue
                sequence = line.split("\t")[9]
                tag = category(sequence, args.motif)
                if tag:
                    files[tag].write(line)

        for tag, sam in paths.items():
            if sam.stat().st_size:
                sort_bam(sam, args.output / f"{bam.stem}.{tag}.sorted.bam", log)
    return args.output


def motif_positions(sequence, motif):
    position = sequence.find(motif)
    while position != -1:
        yield position
        position = sequence.find(motif, position + 1)


def category(sequence, motif):
    if sequence == "*":
        return None
    if sequence.startswith(motif):
        return "start"
    if sequence.endswith(motif):
        return "end"
    for position in motif_positions(sequence, motif):
        right_length = len(sequence) - position - len(motif)
        if position >= K and right_length >= K:
            return "middle20"
    return None


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


def read_lines(command, log):
    command = list(map(str, command))
    with subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=log, text=True
    ) as process:
        try:
            yield from process.stdout
            code = process.wait()
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
                process.wait()
    if code:
        raise subprocess.CalledProcessError(code, command)


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
