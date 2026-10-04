#!/usr/bin/env python3
"""FASTA фланков и участков между ними."""

import argparse
import csv
import os
import shutil
import subprocess
from contextlib import ExitStack
from pathlib import Path


def between_table(table, out, reference, log):
    _, rows = read_table(table)
    folder = out / table.stem
    folder.mkdir()
    required = {
        "query_id",
        "query_seq",
        "query_len",
        "aln_chr",
        "aln_start0",
        "aln_end0",
    }
    if rows and not required <= rows[0].keys():
        return
    best = write_flanks(rows, folder)
    write_between(best, folder, reference, log)


def write_flanks(rows, folder):
    best = {"L": {}, "R": {}}
    seen = {"L": set(), "R": set()}
    with open(folder / "left20.fa", "w") as left, open(
        folder / "right20.fa", "w"
    ) as right:
        files = {"L": left, "R": right}
        for row in rows:
            query = row.get("query_id", "").strip()
            sequence = row.get("query_seq", "").strip()
            read_id, _, side = read_key(query)
            if (
                number(row.get("query_len"), None) != 20
                or len(sequence) != 20
                or side not in files
            ):
                continue
            if sequence not in seen[side]:
                files[side].write(f">{query}\n{sequence}\n")
                seen[side].add(sequence)
            current = best[side].get(read_id)
            if current is None or coordinate_key(row) < coordinate_key(current):
                best[side][read_id] = row
    return best


def write_between(best, folder, reference, log):
    seen = set()
    files = {}
    paired_reads = sorted(best["L"].keys() & best["R"].keys())
    with ExitStack() as stack:
        for read_id in paired_reads:
            left = best["L"][read_id]
            right = best["R"][read_id]
            if left["aln_chr"] != right["aln_chr"]:
                continue
            start = number(left["aln_end0"], None)
            end = number(right["aln_start0"], None)
            if start is None or end is None or end <= start:
                continue
            length = end - start
            chrom = left["aln_chr"]
            sequence = reference_sequence(reference, chrom, start, end, log)
            if len(sequence) != length or sequence in seen:
                continue
            seen.add(sequence)
            if length not in files:
                files[length] = stack.enter_context(
                    open(folder / f"between_len_{length}.fa", "w")
                )
            files[length].write(
                f">{read_id}|{chrom}:{start}-{end}|len={length}\n{sequence}\n"
            )


def reference_sequence(reference, chrom, start, end, log):
    result = subprocess.run(
        ["samtools", "faidx", reference, f"{chrom}:{start + 1}-{end}"],
        stdout=subprocess.PIPE,
        stderr=log,
        text=True,
    )
    if result.returncode:
        return ""
    return "".join(
        line.strip() for line in result.stdout.splitlines() if not line.startswith(">")
    )


def run(args, log):
    source = args.input
    out = args.output
    reference = out / "reference" / args.reference.name
    reference.parent.mkdir()
    shutil.copy2(args.reference, reference)
    subprocess.run(["samtools", "faidx", reference], stdout=log, stderr=log, check=True)
    directories = args.tables or [source if source.is_dir() else source.parent]
    tables = sorted(
        {
            table
            for folder in directories
            for pattern in ("*.csv", "*.tsv")
            for table in folder.glob(pattern)
        }
    )
    if not tables:
        raise ValueError("No CSV/TSV tables for stage 8")
    if len({p.stem for p in tables}) != len(tables):
        raise ValueError("Tables have duplicate output names")
    for table in tables:
        between_table(table, out, reference, log)
    return out


def number(value, default=10**18):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def read_key(query):
    parts = (query or "").split("|") + ["", ""]
    return parts[0], parts[1], parts[2]


def coordinate_key(row):
    return (
        number(row.get("aln_start0")),
        number(row.get("aln_end0")),
        row.get("aln_chr", ""),
    )


def read_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        if path.suffix.lower() == ".csv":
            reader = csv.DictReader(handle)
            return reader.fieldnames, list(reader)

        lines = (line.rstrip("\n") for line in handle if line.strip())
        fields = next(lines, "").split("\t")
        rows = []
        for line in lines:
            values = line.split("\t")
            values.extend([""] * (len(fields) - len(values)))
            rows.append(dict(zip(fields, values)))
        return fields, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--tables", type=Path, nargs="+")
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
