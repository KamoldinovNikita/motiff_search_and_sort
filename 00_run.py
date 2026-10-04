#!/usr/bin/env python3
"""Последовательный запуск этапов 01–09."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

MOTIF = "GCCAAGGC"

STAGES = {
    1: ("01_trim.py", "01_trimmed_fastp"),
    2: ("02_select_motif.py", f"02_motif_{MOTIF}"),
    3: ("03_align_reads.py", "03_bowtie2_align"),
    4: ("04_remove_regions.py", "04_filtered_no_regions"),
    5: ("05_split_reads.py", f"05_split_{MOTIF}"),
    6: ("06_align_flanks.py", "06_kmers20_align_gene_only"),
    7: ("07_select_candidates.py", "07_middle_candidates"),
    8: ("08_extract_sequences.py", "08_between_fasta"),
    9: ("09_make_logos.py", "09_logos_unique"),
}

HERE = Path(__file__).resolve().parent


PACKAGES = (
    "python=3.11",
    "pip",
    "fastp=1.0.1",
    "seqkit=2.12.0",
    "bowtie2=2.5.5",
    "samtools=1.21",
    "bedtools=2.31.1",
    "ghostscript",
)


def ask_path(label, value=None):
    value = str(value) if value is not None else input(f"{label}: ")
    return Path(value.strip().strip('"').strip("'")).expanduser().resolve()


def install_packages(prefix, log, conda=None):
    manager = str(conda) if conda else os.environ.get("CONDA_EXE") or "conda"
    action = "install" if prefix.exists() else "create"
    print(f"Установка пакетов. Журнал: {log}", flush=True)
    with open(log, "a") as handle:
        subprocess.run(
            [
                manager,
                action,
                "--yes",
                "--prefix",
                str(prefix),
                "--override-channels",
                "--channel",
                "conda-forge",
                "--channel",
                "bioconda",
                "--strict-channel-priority",
                *PACKAGES,
            ],
            stdout=handle,
            stderr=handle,
            check=True,
        )
        subprocess.run(
            [
                str(prefix / "bin/python"),
                "-m",
                "pip",
                "--isolated",
                "install",
                "weblogo==3.9.0",
            ],
            stdout=handle,
            stderr=handle,
            check=True,
        )
    return prefix / "bin"


def run_pipeline(args):
    if args.start_stage > args.stop_stage:
        raise ValueError("Начальный этап не может быть больше последнего")
    stages = range(args.start_stage, args.stop_stage + 1)
    args.input = ask_path("Входной файл или папка", args.input)
    if any(number in stages for number in (3, 6, 8)):
        args.reference = ask_path("Референс FASTA", args.reference)
    if 6 in stages:
        args.annotation = ask_path("Аннотация GFF3", args.annotation)
    args.output = ask_path("Новая папка результатов", args.output)
    args.output.mkdir(parents=True)

    install = args.install_packages == "yes"
    if args.install_packages == "ask":
        answer = input("Установить пакеты в отдельную среду рядом с кодом? [д/Н]: ")
        install = answer.strip().lower() in ("д", "да", "y", "yes")
    environment = HERE / ".environment"
    directories = list(args.bin_dir)
    if install:
        directories.insert(
            0, install_packages(environment, args.output / "install.log", args.conda)
        )
    elif (environment / "conda-meta/history").is_file():
        directories.append(environment / "bin")
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(
        [str(path.resolve()) for path in directories] + [env["PATH"]]
    )

    source = args.input
    completed = {}
    for number in stages:
        script, folder = STAGES[number]
        target = args.output / folder
        command = [
            sys.executable,
            str(HERE / script),
            "--input",
            str(source),
            "--output",
            str(target),
        ]
        if number in (2, 5, 6):
            command += ["--motif", MOTIF]
        if number in (3, 6, 8):
            command += ["--reference", str(args.reference)]
        if number == 6:
            command += ["--annotation", str(args.annotation)]
        if number == 7:
            command += ["--mode", str(args.mode)]
        if number == 8:
            folders = [completed[n] for n in (6, 7) if n in completed]
            if args.start_stage == 7:
                folders.append(args.input if args.input.is_dir() else args.input.parent)
            if folders:
                command += ["--tables", *map(str, folders)]
        if number == 9 and 8 in completed:
            command.append("--recursive")
        print(f"{number}/9 {script}", flush=True)
        subprocess.run(command, env=env, check=True)
        completed[number] = target
        source = target
    print(f"Результаты: {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--annotation", type=Path)
    parser.add_argument(
        "--install-packages", choices=("ask", "yes", "no"), default="ask"
    )
    parser.add_argument("--conda", type=Path)
    parser.add_argument("--bin-dir", action="append", type=Path, default=[])
    parser.add_argument("--mode", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--start-stage", type=int, choices=range(1, 10), default=1)
    parser.add_argument("--stop-stage", type=int, choices=range(1, 10), default=9)
    args = parser.parse_args()
    try:
        run_pipeline(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Ошибка: {error}\n")


if __name__ == "__main__":
    main()
