#!/usr/bin/env python3
"""Построение логотипов."""

import argparse
import os
import subprocess
from pathlib import Path


def fasta_size(path):
    count = 0
    length = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
            elif count == 1:
                length += len(line.rstrip("\n"))
    return count, length


def run(args, log):
    source = args.input
    out = args.output
    search = source.rglob if args.recursive else source.glob
    files = sorted(search("*.fa")) + sorted(search("*.fasta"))
    if args.recursive:
        files = [f for f in files if f.relative_to(source).parts[0] != "reference"]
    if not files:
        raise ValueError("No FASTA files for stage 9")
    for fasta in files:
        folder = out / fasta.parent.relative_to(source)
        temporary = folder / "tmp"
        temporary.mkdir(parents=True, exist_ok=True)
        unique = temporary / f"{fasta.stem}.uniq.fa"
        run_command(["seqkit", "rmdup", "-s", fasta, "-o", unique], log=log)
        count, length = fasta_size(unique)
        if not count:
            continue
        target = folder / f"{fasta.stem}.n_{count}.logo.png"
        env = {**os.environ, "PYTHONWARNINGS": "ignore::UserWarning"}
        run_command(
            [
                "weblogo",
                "-f",
                unique,
                "-D",
                "fasta",
                "-A",
                "dna",
                "-U",
                "bits",
                "-c",
                "classic",
                "-s",
                "large",
                "--resolution",
                "600",
                "-n",
                str(length),
                "-X",
                "YES",
                "-Y",
                "YES",
                "--errorbars",
                "NO",
                "--scale-width",
                "YES",
                "-F",
                "png",
                "-o",
                target,
            ],
            env=env,
            log=log,
        )
        if not target.is_file() or not target.stat().st_size:
            raise RuntimeError(f"WebLogo did not create {target}")
    return out


def run_command(command, log, env=None):
    subprocess.run(command, stdout=log, stderr=log, check=True, env=env)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--recursive", action="store_true")
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
