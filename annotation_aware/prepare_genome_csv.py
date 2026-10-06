#!/usr/bin/env python3
"""Download NCBI genome/annotation files and create an input CSV for `iss generate-annotated`."""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Optional, Tuple
import pandas as pd
import requests
from Bio import Entrez


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Download genomic FASTA and GTF/GFF files for assemblies listed in a CSV, "
	    "then write the path CSV used by `iss generate-annotated`."
        )
    )
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--genomes-dir", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--email", required=True, help="Email required by NCBI Entrez")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--sleep", type=float, default=0.34)
    args = parser.parse_args()

    Entrez.email = args.email

    input_csv = Path(args.input_csv)
    genomes_dir = Path(args.genomes_dir)
    output_csv = Path(args.output_csv)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(input_csv)
    required = {"organism_name", "assembly_accession"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    rows = []
    for index, row in table.iterrows():
        organism = clean(row["organism_name"])
        accession = clean(row["assembly_accession"])
        if not organism:
            raise ValueError(f"Row {index} has an empty organism_name")
        if not accession:
            raise ValueError(f"Row {index} has an empty assembly_accession")
        accession_no_version = strip_version(accession)

        result = {
            "organism_name": organism,
            "assembly_accession": accession,
            "fna_path": "",
            "gtf_path": "",
            "gff_path": "",
            "download_source": "",
            "download_status": "",
            "download_note": "",
        }

        base_url, source, note = resolve_assembly(accession)
        if not base_url:
            result["download_status"] = "resolve_failed"
            result["download_note"] = note
            rows.append(result)
            print(f"[FAIL] {organism}: {note}", flush=True)
            continue

        assembly_name = base_url.rstrip("/").split("/")[-1]
        assembly_dir = genomes_dir / accession_no_version

        fna = assembly_dir / f"{assembly_name}_genomic.fna.gz"
        gtf = assembly_dir / f"{assembly_name}_genomic.gtf.gz"
        gff = assembly_dir / f"{assembly_name}_genomic.gff.gz"

        fna_ok, fna_note = download(f"{base_url}/{assembly_name}_genomic.fna.gz", fna, args.timeout)
        gtf_ok, gtf_note = download(f"{base_url}/{assembly_name}_genomic.gtf.gz", gtf, args.timeout)

        gff_ok = False
        gff_note = "not attempted"
        if not gtf_ok:
            gff_ok, gff_note = download(
                f"{base_url}/{assembly_name}_genomic.gff.gz", gff, args.timeout
            )

        result["download_source"] = source or ""
        result["fna_path"] = str(fna.resolve()) if fna_ok else ""
        result["gtf_path"] = str(gtf.resolve()) if gtf_ok else ""
        result["gff_path"] = str(gff.resolve()) if gff_ok else ""

        if fna_ok and (gtf_ok or gff_ok):
            result["download_status"] = "ok"
        elif not fna_ok:
            result["download_status"] = "missing_fasta"
        else:
            result["download_status"] = "missing_annotation"

        result["download_note"] = (f"FNA: {fna_note}; GTF: {gtf_note}; GFF: {gff_note}")
        rows.append(result)
        print(f"[{result['download_status'].upper()}] {organism} | {accession}", flush=True,)
        time.sleep(args.sleep)

    pd.DataFrame(rows).to_csv(output_csv, index=False)
    print(f"[DONE] wrote {output_csv}", flush=True)


def clean(value) -> str:
    return "" if pd.isna(value) else str(value).strip()


def strip_version(accession: str) -> str:
    accession = clean(accession)
    if accession.startswith("RS_"):
        accession = accession[3:]
    elif accession.startswith("GB_"):
        accession = accession[3:]
    return accession.split(".")[0]


def resolve_assembly(accession: str) -> Tuple[Optional[str], Optional[str], str]:
    """Return the preferred NCBI assembly HTTPS directory and source database."""
    try:
        with Entrez.esearch(
            db="assembly",
            term=f'"{accession}"[Assembly Accession]',
            retmode="xml",
        ) as handle:
            result = Entrez.read(handle)

        ids = result.get("IdList", [])
        if not ids:
            return None, None, f"No NCBI Assembly record found for {accession}"

        with Entrez.esummary(db="assembly", id=ids[0], retmode="xml") as handle:
            summary = Entrez.read(handle)

        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        refseq = clean(doc.get("FtpPath_RefSeq", ""))
        genbank = clean(doc.get("FtpPath_GenBank", ""))

        if refseq:
            return refseq.replace("ftp://", "https://"), "RefSeq", ""
        if genbank:
            return genbank.replace("ftp://", "https://"), "GenBank", ""
        return None, None, f"No downloadable assembly path found for {accession}"
    except Exception as exc:
        return None, None, f"Failed to resolve {accession}: {exc}"


def download(url: str, destination: Path, timeout: int) -> Tuple[bool, str]:
    if destination.is_file() and destination.stat().st_size > 0:
        return True, "already present"
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(url, stream=True, timeout=timeout) as response:
            if response.status_code != 200:
                return False, f"HTTP {response.status_code}"
            with destination.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        output.write(chunk)
        return True, "downloaded"
    except Exception as exc:
        if destination.exists():
            destination.unlink()
        return False, str(exc)


if __name__ == "__main__":
    main()
