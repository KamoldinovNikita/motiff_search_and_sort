#!/usr/bin/env python3
"""Выбор совпадений L и R."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


TABLE_PREFIX = "kmers20_alignments_gene_only"

CHIMERA_TABLE = "tails_alignments_gene_only.grouped_only_chim.tsv"


def run(args):
    table = args.input
    if table.is_dir():
        table = table / f"{TABLE_PREFIX}.grouped.tsv"
    fields, rows = read_table(table)
    groups = defaultdict(lambda: {"L": [], "R": []})
    for row in rows:
        read_id, tag, side = read_key(row.get("query_id"))
        if tag == "middle" and side in ("L", "R"):
            groups[read_id][side].append(row)

    kept = []
    for read_id in sorted(groups):
        left = choose_best(groups[read_id]["L"], args.mode)
        right = choose_best(groups[read_id]["R"], args.mode)
        left_gene = (left.get("gff_gene") or "").strip() if left else ""
        right_gene = (right.get("gff_gene") or "").strip() if right else ""
        if left_gene and left_gene == right_gene:
            continue
        kept.extend(row for row in (left, right) if row is not None)
    write_grouped(args.output / CHIMERA_TABLE, fields, kept)
    return args.output


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


def row_key(row):
    read_id, _, side = read_key(row.get("query_id"))
    return (
        read_id,
        {"L": 0, "R": 1}.get(side, 2),
        row.get("query_id", ""),
        row.get("aln_chr", ""),
        number(row.get("aln_start0")),
    )


def choose_best(rows, mode=2):
    if not rows:
        return None
    if mode == 3:
        return max(
            rows,
            key=lambda row: (
                number(row.get("mapq"), -1),
                -number(row.get("aln_start0")),
                row.get("aln_chr", ""),
            ),
        )
    if mode == 2:
        gene_hits = [row for row in rows if (row.get("gff_gene") or "").strip()]
        if gene_hits:
            rows = gene_hits
    return min(rows, key=coordinate_key)


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


def write_grouped(path, fields, rows):
    rows.sort(key=row_key)
    previous = None
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            read_id, _, side = read_key(row.get("query_id"))
            if previous is not None:
                if previous[0] != read_id:
                    handle.write("\n\n")
                elif previous[1] != side:
                    handle.write("\n")
            handle.write("\t".join(row.get(field, "") for field in fields) + "\n")
            previous = read_id, side


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", type=int, choices=(1, 2, 3), default=2)
    args = parser.parse_args()
    args.input = args.input.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True)
    run(args)


if __name__ == "__main__":
    main()
