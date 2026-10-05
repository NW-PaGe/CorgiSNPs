#!/usr/bin/env python3
"""
Convert a filtered VCF to a CorgiSNPs-like genome-wide wide CSV.

Python: 3.8+

Output format:
    CHROM,POS,FILTER,<reference-column>,<sample1>,<sample2>,...

One row is written for every position in the supplied reference FASTA/FNA.
Each sample column contains:

- the alternate allele for passing SNP calls;
- N for ambiguous, filtered, or below-minimum-depth calls;
- the reference allele otherwise.

FILTER values produced by this script:

- constant: position was not represented by an included SNP-like VCF record;
- ambiguous_reference: reference FASTA contains N at this position;
- above_max_amb_samples: number of N calls at this position exceeds
  --max_amb_samples.

Multiple FILTER values are separated with semicolons.
"""

import argparse
from collections import defaultdict
import csv
import sys

import vcfTools


DEFAULT_MAX_AMBIGUOUS = 1_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert VCF SNP calls to a CorgiSNPs-like genome-wide wide CSV."
    )
    parser.add_argument("infile", help="VCF file")
    parser.add_argument(
        "--reference",
        "--reference-fasta",
        required=True,
        dest="reference",
        help="Reference FASTA/FNA used for coordinate order and constant/reference calls",
    )
    parser.add_argument(
        "--max_amb_samples",
        type=int,
        default=None,
        help="Maximum number of N calls before marking FILTER=above_max_amb_samples",
    )
    parser.add_argument(
        "--min_depth",
        type=int,
        default=10,
        help='Replace SNP with "N" if depth is less than minimum (default: 10)',
    )
    parser.add_argument(
        "--reference-column-name",
        default="REF",
        help='Name to use for the reference-base column. Default: "REF"',
    )
    parser.add_argument(
        "--delimiter",
        default=",",
        help='Output delimiter. Default: comma. Use "\\t" for tab-delimited output.',
    )
    return parser.parse_args()


def normalize_delimiter(value):
    if value == r"\t":
        return "\t"
    return value


def read_reference_fasta(path):
    """Return ordered contig names and uppercase sequence strings from FASTA/FNA."""
    contig_order = []
    sequences = {}
    current_name = None
    current_parts = []

    with open(path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_parts).upper()

                current_name = line[1:].strip().split()[0]
                if not current_name:
                    raise ValueError(
                        "{}:{} has a FASTA header with no contig name".format(
                            path, line_number
                        )
                    )
                if current_name in sequences or current_name in contig_order:
                    raise ValueError(
                        "{}:{} duplicate FASTA contig name: {}".format(
                            path, line_number, current_name
                        )
                    )
                contig_order.append(current_name)
                current_parts = []
            else:
                if current_name is None:
                    raise ValueError(
                        "{}:{} contains sequence before the first FASTA header".format(
                            path, line_number
                        )
                    )
                current_parts.append(line.strip())

    if current_name is not None:
        sequences[current_name] = "".join(current_parts).upper()

    if not contig_order:
        raise ValueError("{} appears to contain no FASTA records".format(path))

    return contig_order, sequences


def sample_names_for_vcf(header, vcf_path):
    """Return output sample names and names used to look up columns in the VCF."""
    samples = header.get_samples()

    # Kept for compatibility with older vcfTools behavior. Current VcfHeader already
    # replaces ['SAMPLE'] with the VCF path and maps it back to column 9.
    if samples == ["SAMPLE"]:
        sys.stderr.write("No sample name in {}, using file name.\n".format(vcf_path))
        return [vcf_path], {vcf_path: "SAMPLE"}

    return samples, {sample: sample for sample in samples}


def genotype_depth(record, sample_index):
    format_keys = record.format.split(":")
    sample_values = record.genotypes[sample_index].split(":")
    return int(dict(zip(format_keys, sample_values))["DP"])


def is_snp_like_record(record):
    """Return True when the record is a single-base REF with single-base ALT(s)."""
    if len(record.get_ref()) != 1:
        return False
    return any(len(alt) == 1 for alt in record.get_alt_field().split(","))


def record_call(record, caller, genotype, sample_index, min_depth):
    """Return (call, include_site) for this sample.

    call is the base to store for this sample, or None for reference/no-variant.
    include_site indicates whether this sample makes the position eligible as a
    SNP-like VCF site, even if the stored call is N.
    """
    if not record.is_passing(caller):
        return "N", False

    variant_type = record.get_variant_type(caller, genotype)
    if not variant_type:
        return None, False

    if variant_type == "SNP":
        if genotype_depth(record, sample_index) >= min_depth:
            return record.get_alt(genotype), True
        return "N", True

    if variant_type == "uncalled_ambiguous" and is_snp_like_record(record):
        return "N", True

    return "N", False


def collect_calls(vcf_path, min_depth):
    header = vcfTools.VcfHeader(vcf_path)
    caller = header.get_caller()
    samples, sample_lookup = sample_names_for_vcf(header, vcf_path)

    vcf_refs = defaultdict(dict)                     # vcf_refs[chrom][pos] = ref_base in VCF
    calls = defaultdict(lambda: defaultdict(dict))    # calls[sample][chrom][pos] = base or N
    included_positions = defaultdict(set)             # SNP-like positions from VCF
    ambiguous_counts = defaultdict(lambda: defaultdict(int))

    with open(vcf_path, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue

            record = vcfTools.VcfRecord(line)
            chrom = record.get_chrom()
            pos = int(record.get_pos())

            for sample in samples:
                lookup_name = sample_lookup[sample]
                sample_index = header.get_sample_index(lookup_name)
                genotype = record.get_genotype(index=sample_index, min_gq=0)
                call, include_site = record_call(
                    record, caller, genotype, sample_index, min_depth
                )

                if call is None:
                    continue

                vcf_refs[chrom][pos] = record.get_ref().upper()
                calls[sample][chrom][pos] = call.upper() if isinstance(call, str) else call
                if include_site:
                    included_positions[chrom].add(pos)
                if call == "N":
                    ambiguous_counts[chrom][pos] += 1

    return samples, vcf_refs, calls, included_positions, ambiguous_counts


def warn_for_vcf_reference_mismatches(vcf_refs, reference_sequences):
    missing_from_reference = 0
    mismatched_reference = 0

    for chrom in sorted(vcf_refs):
        sequence = reference_sequences.get(chrom)
        if sequence is None:
            missing_from_reference += len(vcf_refs[chrom])
            continue

        for pos, vcf_ref in vcf_refs[chrom].items():
            if pos < 1 or pos > len(sequence):
                missing_from_reference += 1
                continue
            fasta_ref = sequence[pos - 1].upper()
            if fasta_ref != vcf_ref:
                mismatched_reference += 1

    if missing_from_reference:
        sys.stderr.write(
            "WARNING: {} VCF position(s) are not present in the reference FASTA and will not be emitted.\n".format(
                missing_from_reference
            )
        )
    if mismatched_reference:
        sys.stderr.write(
            "WARNING: {} VCF position(s) have REF bases that differ from the reference FASTA. FASTA bases are used in the output REF column.\n".format(
                mismatched_reference
            )
        )


def build_filter_value(fasta_ref, has_included_vcf_site, ambiguous_count, max_ambiguous):
    filters = []

    if ambiguous_count > max_ambiguous:
        filters.append("above_max_amb_samples")
    elif fasta_ref == "N":
        filters.append("ambiguous_reference")
    elif not has_included_vcf_site:
        filters.append("constant")

    return ";".join(filters)


def iter_wide_rows(
    samples,
    reference_contigs,
    reference_sequences,
    calls,
    included_positions,
    ambiguous_counts,
    max_ambiguous,
):
    for chrom in reference_contigs:
        sequence = reference_sequences[chrom]
        for offset, fasta_ref in enumerate(sequence):
            pos = offset + 1
            has_included_vcf_site = pos in included_positions.get(chrom, set())
            filter_value = build_filter_value(
                fasta_ref=fasta_ref,
                has_included_vcf_site=has_included_vcf_site,
                ambiguous_count=ambiguous_counts[chrom][pos],
                max_ambiguous=max_ambiguous,
            )

            row = [chrom, pos, filter_value, fasta_ref]
            for sample in samples:
                row.append(calls[sample][chrom].get(pos, fasta_ref))
            yield row


def main():
    args = parse_args()
    vcf_path = args.infile.rstrip()
    sys.stderr.write("Searching {}\n".format(vcf_path))
    sys.stderr.write("Reading reference {}\n".format(args.reference))

    max_ambiguous = args.max_amb_samples
    if max_ambiguous is None:
        max_ambiguous = DEFAULT_MAX_AMBIGUOUS

    reference_contigs, reference_sequences = read_reference_fasta(args.reference)
    samples, vcf_refs, calls, included_positions, ambiguous_counts = collect_calls(
        vcf_path=vcf_path,
        min_depth=args.min_depth,
    )
    warn_for_vcf_reference_mismatches(vcf_refs, reference_sequences)

    delimiter = normalize_delimiter(args.delimiter)
    writer = csv.writer(sys.stdout, delimiter=delimiter, lineterminator="\n")
    writer.writerow(["CHROM", "POS", "FILTER", args.reference_column_name] + samples)
    writer.writerows(
        iter_wide_rows(
            samples=samples,
            reference_contigs=reference_contigs,
            reference_sequences=reference_sequences,
            calls=calls,
            included_positions=included_positions,
            ambiguous_counts=ambiguous_counts,
            max_ambiguous=max_ambiguous,
        )
    )


if __name__ == "__main__":
    main()
