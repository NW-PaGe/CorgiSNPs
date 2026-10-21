#!/usr/bin/env python3
"""
Run MycoSNP VCF -> wide CSV conversion, then compare against CorgiSNPs full.csv.

Python: 3.8+

Accepts gzipped or plain reference FASTA and MycoSNP VCF. Sets
--max_amb_samples on the converter to floor(10% of VCF sample count),
minimum 1. Reads the CorgiSNPs reference-column name from full.csv and
passes it to the converter as --reference-column-name.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Optional, TextIO


METADATA_COLUMNS = {"CHROM", "POS", "FILTER"}
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
VCF_TO_CSV = os.path.join(SCRIPT_DIR, "vcfSnpsToFasta_genomewide_fullcsv.py")
COMPARE = os.path.join(SCRIPT_DIR, "compare_mycosnp_corgisnps.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a MycoSNP VCF to a genome-wide wide CSV, then compare "
            "against a CorgiSNPs full.csv."
        )
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA/FNA (plain or .gz).",
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="MycoSNP VCF (plain or .gz).",
    )
    parser.add_argument(
        "--corgisnps",
        required=True,
        help="CorgiSNPs full.csv.",
    )
    parser.add_argument(
        "--mycosnp-csv",
        default="mycosnp_full.csv",
        help='Path for the intermediate MycoSNP wide CSV. Default: "mycosnp_full.csv".',
    )
    parser.add_argument(
        "--output",
        default="snp_comparison_summary.csv",
        help='Comparison summary CSV. Default: "snp_comparison_summary.csv".',
    )
    parser.add_argument(
        "--details-output",
        default="snp_comparison_details.csv",
        help='Per-difference details CSV. Default: "snp_comparison_details.csv".',
    )
    parser.add_argument(
        "--min_depth",
        type=int,
        default=10,
        help="Minimum depth for SNP calls in the converter (default: 10).",
    )
    return parser.parse_args()


def is_gzip_path(path: str) -> bool:
    if path.endswith(".gz"):
        return True
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path: str) -> TextIO:
    if is_gzip_path(path):
        return gzip.open(path, "rt")
    return open(path, "r")


def decompress_if_needed(path: str, dest_dir: str, label: str) -> str:
    """Return a plain-text path; decompress into dest_dir when input is gzipped."""
    if not is_gzip_path(path):
        return path

    base = os.path.basename(path)
    if base.endswith(".gz"):
        base = base[:-3]
    out_path = os.path.join(dest_dir, "{}_{}".format(label, base))
    sys.stderr.write("Decompressing {} -> {}\n".format(path, out_path))
    with gzip.open(path, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def count_vcf_samples(vcf_path: str) -> int:
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    raise SystemExit(
                        "VCF header has no sample columns: {}".format(vcf_path)
                    )
                return len(fields) - 9
            if not line.startswith("#"):
                break
    raise SystemExit("Could not find #CHROM header line in {}".format(vcf_path))


def max_amb_samples_from_count(sample_count: int) -> int:
    if sample_count < 1:
        raise SystemExit("VCF has zero samples")
    return max(1, math.floor(sample_count * 0.1))


def read_corgisnps_ref_column_name(path: str) -> str:
    """Return the first non-metadata column after CHROM/POS/FILTER."""
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit("CorgiSNPs file is empty: {}".format(path))

    missing = [column for column in ("CHROM", "POS", "FILTER") if column not in header]
    if missing:
        raise SystemExit(
            "{} is missing required column(s): {}".format(path, ", ".join(missing))
        )

    non_metadata = [column for column in header if column not in METADATA_COLUMNS]
    if not non_metadata:
        raise SystemExit(
            "{} has no reference/sample columns after CHROM, POS, FILTER".format(path)
        )
    return non_metadata[0]


def run_checked(cmd: list, stdout: Optional[TextIO] = None) -> None:
    sys.stderr.write("Running: {}\n".format(" ".join(cmd)))
    subprocess.run(cmd, check=True, stdout=stdout)


def main() -> None:
    args = parse_args()

    for path, label in (
        (args.reference, "reference"),
        (args.vcf, "vcf"),
        (args.corgisnps, "corgisnps"),
    ):
        if not os.path.isfile(path):
            raise SystemExit("{} not found: {}".format(label, path))

    sample_count = count_vcf_samples(args.vcf)
    max_amb = max_amb_samples_from_count(sample_count)
    ref_column_name = read_corgisnps_ref_column_name(args.corgisnps)
    sys.stderr.write(
        "VCF samples: {}; max_amb_samples: {} (10% floored, min 1)\n".format(
            sample_count, max_amb
        )
    )
    sys.stderr.write(
        "CorgiSNPs reference column: {}\n".format(ref_column_name)
    )

    with tempfile.TemporaryDirectory(prefix="mycosnp_corgi_") as tmpdir:
        reference = decompress_if_needed(args.reference, tmpdir, "reference")
        vcf = decompress_if_needed(args.vcf, tmpdir, "vcf")

        convert_cmd = [
            sys.executable,
            VCF_TO_CSV,
            vcf,
            "--reference",
            reference,
            "--max_amb_samples",
            str(max_amb),
            "--min_depth",
            str(args.min_depth),
            "--reference-column-name",
            ref_column_name,
        ]
        mycosnp_csv_dir = os.path.dirname(os.path.abspath(args.mycosnp_csv))
        if mycosnp_csv_dir:
            os.makedirs(mycosnp_csv_dir, exist_ok=True)

        with open(args.mycosnp_csv, "w") as out_handle:
            run_checked(convert_cmd, stdout=out_handle)
        sys.stderr.write("Wrote {}\n".format(args.mycosnp_csv))

        compare_cmd = [
            sys.executable,
            COMPARE,
            "--mycosnp",
            args.mycosnp_csv,
            "--corgisnps",
            args.corgisnps,
            "--output",
            args.output,
            "--details-output",
            args.details_output,
        ]

        run_checked(compare_cmd)


if __name__ == "__main__":
    main()
