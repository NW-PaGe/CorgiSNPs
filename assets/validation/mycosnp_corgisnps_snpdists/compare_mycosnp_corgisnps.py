#!/usr/bin/env python3
"""
Compare sample-wise SNP calls between MycoSNP and CorgiSNPs wide full.csv-style outputs.

Python: 3.8+

Expected input format for both files:
    CHROM,POS,FILTER,<reference-column>,<sample1>,<sample2>,...

Assumptions:
- Both files are in the same coordinate order.
- The first non-metadata column after CHROM/POS/FILTER is the reference column.
- Sample columns must contain the same sample names between files; order may differ.
- MycoSNP sample order is used in the summary output.

Rules:
- Reference column names must be identical between files.
- Coordinates and reference bases must agree between files.
- Rows with ambiguous_reference in either FILTER are excluded from comparison.
- Rows with above_max_amb_samples in MycoSNP FILTER are treated as N for all MycoSNP samples.
- Rows with below_min_cf in CorgiSNPs FILTER are treated as N for all CorgiSNPs samples.
- A call is an alternate allele when it is not N and differs from the reference base.
- A difference is counted when, for the same sample and coordinate:
    1. MycoSNP has an alternate allele and CorgiSNPs has reference.
    2. CorgiSNPs has an alternate allele and MycoSNP has reference.
    3. Both have alternate alleles, but the alternate bases differ.
- Sample-positions where either call is N are skipped.
"""

import argparse
import csv
import itertools
import logging
from typing import Dict, Iterable, List, Optional, Set, Tuple

METADATA_COLUMNS = {"CHROM", "POS", "FILTER"}
MISSING_CALLS = {"", ".", "N"}
AMBIGUOUS_REFERENCE = "ambiguous_reference"
MYCO_FORCE_N_FILTER = "above_max_amb_samples"
CORGI_FORCE_N_FILTER = "below_min_cf"

SummaryMap = Dict[str, Dict[str, int]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Optimized stream comparison of wide MycoSNP and CorgiSNPs full.csv-style SNP calls."
    )
    parser.add_argument(
        "--mycosnp",
        required=True,
        help="Wide MycoSNP CSV. Expected columns: CHROM, POS, FILTER, reference column, sample columns.",
    )
    parser.add_argument(
        "--corgisnps",
        required=True,
        help="CorgiSNPs full.csv. Expected columns: CHROM, POS, FILTER, reference column, sample columns.",
    )
    parser.add_argument(
        "--output",
        default="snp_comparison_summary.csv",
        help="Output summary CSV path. Default: snp_comparison_summary.csv",
    )
    parser.add_argument(
        "--details-output",
        default=None,
        help="Optional output CSV with one row per difference. This is streamed to disk, but may be very large.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level. Default: INFO",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000000,
        help="Log progress every N rows. Default: 1000000. Use 0 to disable progress logging.",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    logging.error(message)
    raise SystemExit(1)


def normalize_call(value: Optional[str]) -> str:
    if value is None:
        return ""
    return value.strip().upper()


def normalize_filter(value: Optional[str]) -> str:
    if value is None:
        return ""
    return value.strip()


def split_filter_values(filter_value: Optional[str]) -> Set[str]:
    value = normalize_filter(filter_value)
    if not value:
        return set()
    parts = value.replace(",", ";").replace(" ", ";").split(";")
    return {part.strip() for part in parts if part.strip()}


def has_filter_value(filter_value: Optional[str], target: str) -> bool:
    return target in split_filter_values(filter_value)


def is_missing_normalized(call: str) -> bool:
    return call in MISSING_CALLS


def is_alt_normalized(call: str, ref: str) -> bool:
    return call not in MISSING_CALLS and ref not in MISSING_CALLS and call != ref


def require_header_columns(header: List[str], path: str) -> None:
    missing = [column for column in ("CHROM", "POS", "FILTER") if column not in header]
    if missing:
        raise ValueError("{} is missing required column(s): {}".format(path, ", ".join(missing)))


def get_ref_and_sample_columns(header: List[str], path: str) -> Tuple[str, List[str]]:
    non_metadata = [column for column in header if column not in METADATA_COLUMNS]
    if not non_metadata:
        raise ValueError("{} has no reference/sample columns after CHROM, POS, FILTER".format(path))
    ref_column = non_metadata[0]
    sample_columns = non_metadata[1:]
    if not sample_columns:
        raise ValueError("{} has no sample columns after reference column {}".format(path, ref_column))
    return ref_column, sample_columns


def index_by_column(header: List[str]) -> Dict[str, int]:
    return {column: index for index, column in enumerate(header)}


def validate_and_build_indexes(
    myco_header: List[str],
    corgi_header: List[str],
    myco_path: str,
    corgi_path: str,
) -> Tuple[str, List[str], Dict[str, int], Dict[str, int], List[int], List[int]]:
    require_header_columns(myco_header, myco_path)
    require_header_columns(corgi_header, corgi_path)

    myco_ref_column, myco_samples = get_ref_and_sample_columns(myco_header, myco_path)
    corgi_ref_column, corgi_samples = get_ref_and_sample_columns(corgi_header, corgi_path)

    if myco_ref_column != corgi_ref_column:
        fail(
            "Reference column names differ: MycoSNP={!r}, CorgiSNPs={!r}".format(
                myco_ref_column, corgi_ref_column
            )
        )

    myco_sample_set = set(myco_samples)
    corgi_sample_set = set(corgi_samples)
    if myco_sample_set != corgi_sample_set:
        only_myco = sorted(myco_sample_set - corgi_sample_set)
        only_corgi = sorted(corgi_sample_set - myco_sample_set)
        if only_myco:
            logging.error("Samples only in MycoSNP CSV: %s", ", ".join(only_myco))
        if only_corgi:
            logging.error("Samples only in CorgiSNPs CSV: %s", ", ".join(only_corgi))
        raise SystemExit(1)

    if myco_samples != corgi_samples:
        logging.warning(
            "Sample columns contain the same names but are in a different order; comparing by sample name"
        )

    myco_idx = index_by_column(myco_header)
    corgi_idx = index_by_column(corgi_header)

    myco_sample_indexes = [myco_idx[sample] for sample in myco_samples]
    corgi_sample_indexes = [corgi_idx[sample] for sample in myco_samples]

    return (
        myco_ref_column,
        myco_samples,
        myco_idx,
        corgi_idx,
        myco_sample_indexes,
        corgi_sample_indexes,
    )


def parse_pos(value: str, path: str, row_number: int) -> int:
    try:
        return int(value)
    except ValueError:
        raise ValueError("{}:{} has non-integer POS: {}".format(path, row_number, value))


def classify_difference_from_normalized(ref: str, myco_call: str, corgi_call: str) -> Optional[str]:
    if is_missing_normalized(myco_call) or is_missing_normalized(corgi_call):
        return None

    myco_alt = is_alt_normalized(myco_call, ref)
    corgi_alt = is_alt_normalized(corgi_call, ref)

    if myco_alt and not corgi_alt:
        return "mycosnp_alt_only"
    if corgi_alt and not myco_alt:
        return "corgisnps_alt_only"
    if myco_alt and corgi_alt and myco_call != corgi_call:
        return "different_alt_alleles"
    return None


def init_summary(samples: Iterable[str]) -> SummaryMap:
    return {
        sample: {
            "positions_compared": 0,
            "positions_skipped_N": 0,
            "positions_excluded_ambiguous_reference": 0,
            "mycosnp_alt_only": 0,
            "corgisnps_alt_only": 0,
            "different_alt_alleles": 0,
            "total_differences": 0,
        }
        for sample in samples
    }


def increment_all(summary: SummaryMap, key: str) -> None:
    for counts in summary.values():
        counts[key] += 1


def open_details_writer(path: Optional[str]):
    if path is None:
        return None, None

    handle = open(path, "w", newline="")
    fieldnames = [
        "sample",
        "chrom",
        "pos",
        "difference_type",
        "mycosnp_filter",
        "corgisnps_filter",
        "mycosnp_ref",
        "mycosnp_call",
        "corgisnps_ref",
        "corgisnps_call",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    return handle, writer


def format_count_percent(count: int, denominator: int) -> str:
    if denominator <= 0:
        return "{} (NA)".format(count)
    percent = 100.0 * count / denominator
    return "{} ({:.3f}%)".format(count, percent)


def format_summary_row(sample: str, counts: Dict[str, int], genome_positions: int) -> Dict[str, object]:
    return {
        "sample": sample,
        "positions_compared": format_count_percent(
            counts["positions_compared"], genome_positions
        ),
        "positions_skipped_N": format_count_percent(
            counts["positions_skipped_N"], genome_positions
        ),
        "positions_excluded_ambiguous_reference": format_count_percent(
            counts["positions_excluded_ambiguous_reference"], genome_positions
        ),
        "mycosnp_alt_only": counts["mycosnp_alt_only"],
        "corgisnps_alt_only": counts["corgisnps_alt_only"],
        "different_alt_alleles": counts["different_alt_alleles"],
        "total_differences": counts["total_differences"],
    }


def write_summary(path: str, summary: SummaryMap, genome_positions: int) -> None:
    fieldnames = [
        "sample",
        "positions_compared",
        "positions_skipped_N",
        "positions_excluded_ambiguous_reference",
        "mycosnp_alt_only",
        "corgisnps_alt_only",
        "different_alt_alleles",
        "total_differences",
    ]

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample in sorted(summary):
            writer.writerow(format_summary_row(sample, summary[sample], genome_positions))


def validate_row_lengths(
    myco_row: List[str],
    corgi_row: List[str],
    myco_header_len: int,
    corgi_header_len: int,
    row_number: int,
) -> None:
    if len(myco_row) != myco_header_len:
        fail(
            "MycoSNP row {} has {} columns but header has {}".format(
                row_number, len(myco_row), myco_header_len
            )
        )
    if len(corgi_row) != corgi_header_len:
        fail(
            "CorgiSNPs row {} has {} columns but header has {}".format(
                row_number, len(corgi_row), corgi_header_len
            )
        )


def stream_compare(args: argparse.Namespace) -> None:
    with open(args.mycosnp, newline="") as myco_handle, open(args.corgisnps, newline="") as corgi_handle:
        myco_reader = csv.reader(myco_handle)
        corgi_reader = csv.reader(corgi_handle)

        try:
            myco_header = next(myco_reader)
        except StopIteration:
            fail("{} appears to be empty".format(args.mycosnp))
        try:
            corgi_header = next(corgi_reader)
        except StopIteration:
            fail("{} appears to be empty".format(args.corgisnps))

        (
            ref_column,
            samples,
            myco_idx,
            corgi_idx,
            myco_sample_indexes,
            corgi_sample_indexes,
        ) = validate_and_build_indexes(myco_header, corgi_header, args.mycosnp, args.corgisnps)

        myco_chrom_i = myco_idx["CHROM"]
        myco_pos_i = myco_idx["POS"]
        myco_filter_i = myco_idx["FILTER"]
        myco_ref_i = myco_idx[ref_column]

        corgi_chrom_i = corgi_idx["CHROM"]
        corgi_pos_i = corgi_idx["POS"]
        corgi_filter_i = corgi_idx["FILTER"]
        corgi_ref_i = corgi_idx[ref_column]

        summary = init_summary(samples)
        details_handle, details_writer = open_details_writer(args.details_output)

        rows_loaded = 0
        myco_force_n_rows = 0
        corgi_force_n_rows = 0
        ambiguous_reference_rows = 0
        sample_count = len(samples)

        try:
            sentinel = object()
            paired_rows = itertools.zip_longest(myco_reader, corgi_reader, fillvalue=sentinel)
            for row_number, (myco_row, corgi_row) in enumerate(paired_rows, start=2):
                if myco_row is sentinel:
                    fail("MycoSNP CSV ended before CorgiSNPs CSV at data row {}".format(row_number))
                if corgi_row is sentinel:
                    fail("CorgiSNPs CSV ended before MycoSNP CSV at data row {}".format(row_number))

                validate_row_lengths(
                    myco_row=myco_row,
                    corgi_row=corgi_row,
                    myco_header_len=len(myco_header),
                    corgi_header_len=len(corgi_header),
                    row_number=row_number,
                )

                myco_chrom = myco_row[myco_chrom_i].strip()
                corgi_chrom = corgi_row[corgi_chrom_i].strip()
                myco_pos = parse_pos(myco_row[myco_pos_i], args.mycosnp, row_number)
                corgi_pos = parse_pos(corgi_row[corgi_pos_i], args.corgisnps, row_number)

                if myco_chrom != corgi_chrom or myco_pos != corgi_pos:
                    fail(
                        "Coordinate mismatch at data row {}: MycoSNP={}:{} CorgiSNPs={}:{}".format(
                            row_number, myco_chrom, myco_pos, corgi_chrom, corgi_pos
                        )
                    )

                ref = normalize_call(myco_row[myco_ref_i])
                corgi_ref = normalize_call(corgi_row[corgi_ref_i])
                if ref != corgi_ref:
                    fail(
                        "Reference mismatch at {}:{}: MycoSNP={} CorgiSNPs={}".format(
                            myco_chrom, myco_pos, ref, corgi_ref
                        )
                    )

                rows_loaded += 1

                if args.progress_every and rows_loaded % args.progress_every == 0:
                    logging.info("Processed %d rows", rows_loaded)

                myco_filter = normalize_filter(myco_row[myco_filter_i])
                corgi_filter = normalize_filter(corgi_row[corgi_filter_i])

                myco_filter_values = split_filter_values(myco_filter)
                corgi_filter_values = split_filter_values(corgi_filter)

                if MYCO_FORCE_N_FILTER in myco_filter_values:
                    myco_force_n_rows += 1
                if CORGI_FORCE_N_FILTER in corgi_filter_values:
                    corgi_force_n_rows += 1

                if AMBIGUOUS_REFERENCE in myco_filter_values or AMBIGUOUS_REFERENCE in corgi_filter_values:
                    ambiguous_reference_rows += 1
                    increment_all(summary, "positions_excluded_ambiguous_reference")
                    continue

                myco_force_n = MYCO_FORCE_N_FILTER in myco_filter_values
                corgi_force_n = CORGI_FORCE_N_FILTER in corgi_filter_values

                # Fast path: when an entire row is forced to N in either file, every
                # sample-position is skipped and no per-sample allele comparison is needed.
                if myco_force_n or corgi_force_n:
                    increment_all(summary, "positions_skipped_N")
                    continue

                # Fast path: both files mark the position as constant, so every sample
                # matches the reference and no per-sample inspection is needed.
                if (
                    "constant" in myco_filter_values
                    and "constant" in corgi_filter_values
                ):
                    increment_all(summary, "positions_compared")
                    continue

                for i in range(sample_count):
                    sample = samples[i]
                    myco_call = normalize_call(myco_row[myco_sample_indexes[i]])
                    corgi_call = normalize_call(corgi_row[corgi_sample_indexes[i]])

                    counts = summary[sample]

                    if myco_call in MISSING_CALLS or corgi_call in MISSING_CALLS:
                        counts["positions_skipped_N"] += 1
                        continue

                    counts["positions_compared"] += 1

                    difference_type = classify_difference_from_normalized(
                        ref,
                        myco_call,
                        corgi_call,
                    )

                    if difference_type is None:
                        continue

                    counts[difference_type] += 1
                    counts["total_differences"] += 1

                    if details_writer is not None:
                        details_writer.writerow(
                            {
                                "sample": sample,
                                "chrom": myco_chrom,
                                "pos": myco_pos,
                                "difference_type": difference_type,
                                "mycosnp_filter": myco_filter,
                                "corgisnps_filter": corgi_filter,
                                "mycosnp_ref": ref,
                                "mycosnp_call": myco_call,
                                "corgisnps_ref": ref,
                                "corgisnps_call": corgi_call,
                            }
                        )

        finally:
            if details_handle is not None:
                details_handle.close()

    logging.info("Loaded/processed %d paired rows", rows_loaded)
    logging.info("Reference column: %s", ref_column)
    logging.info("Shared samples: %d", len(samples))
    logging.info("MycoSNP rows forced to N by %s: %d", MYCO_FORCE_N_FILTER, myco_force_n_rows)
    logging.info("CorgiSNPs rows forced to N by %s: %d", CORGI_FORCE_N_FILTER, corgi_force_n_rows)
    logging.info("Rows excluded by ambiguous_reference in either file: %d", ambiguous_reference_rows)

    write_summary(args.output, summary, rows_loaded)
    logging.info("Wrote summary: %s", args.output)
    if args.details_output:
        logging.info("Wrote details: %s", args.details_output)


def main() -> None:
    args = parse_args()
    logging.basicConfig(format="%(levelname)s: %(message)s", level=getattr(logging, args.log_level))

    try:
        stream_compare(args)
    except ValueError as error:
        logging.error(str(error))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
