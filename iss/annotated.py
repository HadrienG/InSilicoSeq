#!/usr/bin/env python3

from __future__ import annotations

import gzip
import os
import re
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from iss import generator


def generate_annotated_reads(args):
    """Implementation of `iss generate-annotated`."""

    if args.read_pairs <= 0:
        raise ValueError("--read-pairs must be positive")

    if args.output_read_length is not None and args.output_read_length <= 0:
        raise ValueError("--output-read-length must be positive")

    error_model = generator.load_error_model(
        args.mode,
        args.seed,
        args.model,
        args.fragment_length,
        args.fragment_length_sd,
        True,
    )

    table = pd.read_csv(args.csv)

    required = {"organism_name", "fna_path"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    status_rows = []

    for index, row in table.iterrows():
        organism = clean(row["organism_name"])
        fna_path = clean(row["fna_path"])
        annotation_path = choose_annotation(row)

        if not organism:
            raise ValueError(f"Row {index} has an empty organism_name")

        if not fna_path or not Path(fna_path).is_file():
            raise FileNotFoundError(f"Missing FASTA for {organism}: {fna_path}")

        if (args.reading_frame is not None or args.only_full_coding) and not annotation_path:
            raise FileNotFoundError(
                f"{organism} requires an existing gtf_path "
                f"or gff_path for reading-frame/CDS filtering"
            )

        prefix = safe_name(organism) + frame_tag(args.reading_frame)

        r1 = output_dir / f"{prefix}_R1.fastq"
        r2 = output_dir / f"{prefix}_R2.fastq"
        mutations = output_dir / f"{prefix}_mutations.tsv"
        annotations = (output_dir / f"{prefix}_R1_annotations.tsv")

        final_paths = [r1, r2, mutations, annotations]

        if (all(path.is_file() for path in final_paths) and not args.force):
            print(f"[SKIP] {organism}: outputs already exist", flush=True)

            status_rows.append(
                {
                    "organism_name": organism,
                    "requested_read_pairs": args.read_pairs,
                    "generated_read_pairs": "",
                    "status": "skipped_existing",
                }
            )
            continue

        temporary_outputs = [Path(str(path) + ".tmp") for path in final_paths]

        for path in temporary_outputs:
            if path.exists():
                path.unlink()

        with tempfile.TemporaryDirectory(prefix=f"iss_annotated_{prefix}_") as temp_dir_name:

            temp_dir = Path(temp_dir_name)

            fasta = decompress_if_needed(fna_path, temp_dir, "genome.fasta")

            annotation = ""

            if annotation_path:
                annotation = decompress_if_needed(annotation_path, temp_dir,
                    decompressed_annotation_name(annotation_path))

            try:
                with (
                    temporary_outputs[0].open("w") as r1_handle,
                    temporary_outputs[1].open("w") as r2_handle,
                    temporary_outputs[2].open("w") as mutation_handle,
                    temporary_outputs[3].open("w") as annotation_handle,
                ):

                    annotation_handle.write(
                        "read_id\tchromosome\tstart\tend\t"
                        "is_cds\tfully_contained_in_cds\t"
                        "overlap_pct\tcoding_bases\t"
                        "read_length\treading_frame\t"
                        "cds_strand\tcds_phase\t"
                        "best_cds_id\n"
                    )

                    accepted = (
                        generator.simulate_and_add_info_to_reads(
                            annotation,
                            fasta,
                            error_model,
                            n_pairs=args.read_pairs,
                            cpu_number=index + 1,
                            forward_handle=r1_handle,
                            reverse_handle=r2_handle,
                            mutations_handle=mutation_handle,
                            sequence_type="metagenomics",
                            gc_bias=args.gc_bias,
                            target_reading_frame=(args.reading_frame),
                            only_full_coding=(args.only_full_coding),
                            annotation_handle=(annotation_handle),
                            output_read_length=(args.output_read_length),
                        )
                    )

                accepted = int(accepted or 0)

                if accepted <= 0:
                    raise RuntimeError(
                        f"No read pairs were generated for {organism}"
                    )

                for temporary, final in zip(temporary_outputs, final_paths):
                    os.replace(temporary, final)

                status = ("ok" if accepted == args.read_pairs else "ok_partial")

                status_rows.append(
                    {
                        "organism_name": organism, "requested_read_pairs": (args.read_pairs),
                        "generated_read_pairs": accepted, "status": status,
                    }
                )

                print(
                    f"[{status.upper()}] "
                    f"{organism}: "
                    f"{accepted}/{args.read_pairs} pairs",
                    flush=True,
                )

            except Exception:
                for path in temporary_outputs:
                    if path.exists():
                        path.unlink()
                raise

    status_path = (output_dir / "generation_status.csv")

    pd.DataFrame(status_rows).to_csv(status_path, index=False,)

    print(f"[DONE] wrote {status_path}", flush=True,)


def clean(value):
    return (
        ""
        if pd.isna(value)
        else str(value).strip()
    )


def safe_name(value):
    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        clean(value),
    ).strip("_")


def choose_annotation(row):
    for column in ("gtf_path", "gff_path"):
        path = clean(row.get(column, ""))

        if path and Path(path).is_file():
            return path

    return ""


def decompressed_annotation_name(source):
    name = Path(source).name.lower()

    if ".gtf" in name:
        return "annotation.gtf"

    if ".gff" in name:
        return "annotation.gff"

    return "annotation.txt"


def decompress_if_needed(
    source,
    output_dir,
    output_name,
):
    if not source.endswith(".gz"):
        return source

    destination = output_dir / output_name

    with (
        gzip.open(source, "rb") as input_handle,
        destination.open("wb") as output_handle,
    ):
        shutil.copyfileobj(input_handle, output_handle,)

    return str(destination)


def frame_tag(frame):
    if frame is None:
        return ""

    sign = "p" if frame > 0 else "m"

    return f"_frame_{sign}{abs(frame)}"
