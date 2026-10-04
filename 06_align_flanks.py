#!/usr/bin/env python3
"""Выравнивание и аннотация фланков."""

import argparse
import csv
import os
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


MOTIF = "GCCAAGGC"

K = 20

FIELDS = [
    "source_file",
    "query_id",
    "query_seq",
    "query_len",
    "read_spaced",
    "aln_chr",
    "aln_start0",
    "aln_end0",
    "aln_strand",
    "mapq",
    "cigar",
    "gff_gene",
    "gene_start0",
    "gene_end0",
    "gene_strand",
]

TABLE_PREFIX = "kmers20_alignments_gene_only"


def gene_bed(args, reference, target):
    with open(reference) as handle:
        contig = next(line[1:].split()[0] for line in handle if line.startswith(">"))
    with open(args.annotation) as source, open(target, "w") as out:
        for line in source:
            if line.startswith("#") or not line.strip():
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = line.rstrip(
                "\n"
            ).split("\t")
            if feature != "gene":
                continue
            if chrom in ("NC_000913.3", "chrNC_000913.3", "NC_000913", "chrNC_000913"):
                chrom = contig
            attributes = dict(
                item.split("=", 1)
                for item in attributes.strip().split(";")
                if "=" in item
            )
            gene = attributes.get("gene", "") or attributes.get("gene_name", "")
            locus = attributes.get("locus_tag", "")
            name = attributes.get("Name", "")
            identifier = attributes.get("ID", "")
            label = gene or locus or name or identifier
            fields = [
                chrom,
                str(int(start) - 1),
                end,
                label,
                "0",
                strand,
                label,
                gene,
                locus,
                name,
                identifier,
            ]
            out.write("\t".join(fields) + "\n")


def extract_flanks(bam, tag, target, log, motif):
    motif_length = len(motif)
    with open(target, "w") as out:
        for line in read_lines(["samtools", "view", bam], log):
            fields = line.rstrip("\n").split("\t")
            read_id, sequence = fields[0], fields[9]
            if not sequence or sequence == "*":
                continue

            flanks = []
            if (
                tag == "start"
                and sequence.startswith(motif)
                and len(sequence) >= motif_length + K
            ):
                flank = sequence[motif_length : motif_length + K]
                flanks.append(("R", 1, flank))
            elif (
                tag == "end"
                and sequence.endswith(motif)
                and len(sequence) >= motif_length + K
            ):
                position = len(sequence) - motif_length
                flank = sequence[position - K : position]
                flanks.append(("L", position + 1, flank))
            elif tag == "middle":
                for position in motif_positions(sequence, motif):
                    after_motif = position + motif_length
                    if position >= K and after_motif + K <= len(sequence):
                        flanks.append(
                            ("L", position + 1, sequence[position - K : position])
                        )
                        flanks.append(
                            ("R", position + 1, sequence[after_motif : after_motif + K])
                        )
            for side, position, flank in flanks:
                out.write(f">{read_id}|{tag}|{side}|motifpos{position}|k{K}\n{flank}\n")


def sam_to_bed(source, target):
    with open(source) as sam, open(target, "w") as bed:
        for line in sam:
            if line.startswith("@") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            read_id, flag, chrom, position, mapq, cigar = fields[:6]
            sequence = fields[9]
            if chrom == "*" or int(position) == 0:
                continue
            length = sum(
                int(size)
                for size, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
                if operation in "MDN=X"
            )
            start = int(position) - 1
            strand = "-" if int(flag) & 16 else "+"
            row = [
                chrom,
                str(start),
                str(start + length),
                read_id,
                mapq,
                strand,
                sequence,
                cigar,
            ]
            bed.write("\t".join(row) + "\n")


def annotate(bed, genes, bam, log, motif):
    alignments = [
        line.split("\t") for line in bed.read_text().splitlines() if line.strip()
    ]
    read_ids = {read_key(alignment[3])[0] for alignment in alignments}
    sequences = {}
    if read_ids:
        for line in read_lines(["samtools", "view", bam], log):
            fields = line.rstrip("\n").split("\t")
            read_id, sequence = fields[0], fields[9]
            if read_id in read_ids:
                sequences.setdefault(read_id, sequence)

    def make_row(alignment, gene="", gene_start="", gene_end="", gene_strand=""):
        chrom, start, end, query_id, mapq, strand, sequence, cigar = alignment
        read_id = read_key(query_id)[0]
        return {
            "source_file": bam.name,
            "query_id": query_id,
            "query_seq": sequence,
            "query_len": str(len(sequence)),
            "read_spaced": spaced_read(sequences.get(read_id, ""), motif),
            "aln_chr": chrom,
            "aln_start0": start,
            "aln_end0": end,
            "aln_strand": strand,
            "mapq": mapq,
            "cigar": cigar,
            "gff_gene": gene,
            "gene_start0": gene_start,
            "gene_end0": gene_end,
            "gene_strand": gene_strand,
        }

    annotated = set()
    rows = []
    command = ["bedtools", "intersect", "-wa", "-wb", "-a", bed, "-b", genes]
    for line in read_lines(command, log):
        fields = line.rstrip("\n").split("\t")
        alignment, gene = fields[:8], fields[8:]
        rows.append(make_row(alignment, gene[3], gene[1], gene[2], gene[5]))
        annotated.add(tuple(alignment))
    for alignment in alignments:
        if tuple(alignment) not in annotated:
            rows.append(make_row(alignment))
    return rows


def run(args, log):
    source = args.input
    out = args.output
    temporary = out / "tmp"
    temporary.mkdir()
    reference, index = prepare_reference(args.reference, out, log)
    genes = temporary / "genes.only.bed"
    gene_bed(args, reference, genes)
    rows = []
    for bam in input_files(source, "*.bam"):
        tag = next((t for t in ("start", "middle", "end") if t in bam.name), None)
        if tag is None:
            continue
        base = temporary / f"{bam.name}.{tag}.kmers20"
        fasta, sam, bed = (
            Path(str(base) + suffix) for suffix in (".fa", ".sam", ".bed")
        )
        extract_flanks(bam, tag, fasta, log, args.motif)
        if fasta.stat().st_size:
            run_command(
                [
                    "bowtie2",
                    "-x",
                    index,
                    "-f",
                    "-U",
                    fasta,
                    "-a",
                    "--end-to-end",
                    "-N",
                    "0",
                    "--mp",
                    "999,999",
                    "--rdg",
                    "999,999",
                    "--rfg",
                    "999,999",
                    "-p",
                    "4",
                    "-S",
                    sam,
                ],
                log=log,
            )
        else:
            sam.touch()
        sam_to_bed(sam, bed)
        rows.extend(annotate(bed, genes, bam, log, args.motif))
    write_csv(out / f"{TABLE_PREFIX}.all.csv", FIELDS, rows)
    write_grouped(out / f"{TABLE_PREFIX}.grouped.tsv", FIELDS, rows)
    by_query = defaultdict(list)
    for row in rows:
        by_query[row["query_id"]].append(row)
    best = sorted((choose_best(items) for items in by_query.values()), key=row_key)
    write_csv(out / f"{TABLE_PREFIX}.best_per_query.csv", FIELDS, best)
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


def row_key(row):
    read_id, _, side = read_key(row.get("query_id"))
    return (
        read_id,
        {"L": 0, "R": 1}.get(side, 2),
        row.get("query_id", ""),
        row.get("aln_chr", ""),
        number(row.get("aln_start0")),
    )


def choose_best(rows):
    gene_hits = [row for row in rows if (row.get("gff_gene") or "").strip()]
    return min(gene_hits or rows, key=coordinate_key)


def write_csv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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


def motif_positions(sequence, motif):
    position = sequence.find(motif)
    while position != -1:
        yield position
        position = sequence.find(motif, position + 1)


def spaced_read(sequence, motif):
    if not sequence:
        return "NOT_FOUND"
    start = sequence.find(motif)
    if start == -1:
        return sequence
    end = start + len(motif)
    left_start = max(0, start - K)
    parts = [
        sequence[:left_start],
        sequence[left_start:start],
        motif,
        sequence[end : end + K],
        sequence[end + K :],
    ]
    return " ".join(part for part in parts if part)


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
    parser.add_argument("--motif", default=MOTIF)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
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
