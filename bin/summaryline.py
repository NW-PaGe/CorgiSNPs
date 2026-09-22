#!/usr/bin/env python3
"""
Summarize outputs from bioinformatics workflows.

Aggregates read statistics, species identification, assembly metrics,
and performs automated QC checks.
"""

import json
import csv
import screed
import argparse
import logging
from typing import List, Dict, Any, Optional

# ----------------------------
# Helper Functions
# ----------------------------

def normalize_value(value: Any) -> Optional[str]:
    """Normalize string values for comparison."""
    if value is None:
        return None
    return str(value).lower().strip()


def calculate_z_score(value: float, mean: float, sd: float) -> Optional[float]:
    """Calculate z-score if standard deviation is valid."""
    if sd and sd > 0:
        return (value - mean) / sd
    return None


def compare_values(a: Any, op: str, b: Any) -> bool:
    """Compare two values using the specified operator."""
    operators = {
        "<": lambda x, y: x < y,
        "<=": lambda x, y: x <= y,
        "==": lambda x, y: x == y,
        "!=": lambda x, y: x != y,
        ">=": lambda x, y: x >= y,
        ">": lambda x, y: x > y,
    }
    return operators[op](a, b)


# ----------------------------
# File Loading Functions
# ----------------------------

def load_json(path: str, source: Optional[str] = None) -> Dict[str, Any]:
    """Load JSON file and return as dict."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    logging.info(f"Loaded JSON with {len(data)} top-level keys ({source})")
    return data


def load_csv(path: str, source: Optional[str] = None) -> List[Dict[str, str]]:
    """Load CSV into a list of dictionaries (rows)."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        data = [row for row in reader]
    logging.info(f"Loaded CSV with {len(data)} rows ({source})")
    return data


def load_fasta(path: str):
    """Return a context manager over screed records for streaming."""
    return screed.open(path)


# ----------------------------
# Data Parsing Functions
# ----------------------------

def calculate_genome_stats(records) -> Dict[str, Any]:
    """
    Stream through FASTA and calculate assembly statistics.
    
    Returns:
        Dictionary with contig count, total length, and GC content
    """
    contigs = 0
    total_length = 0
    gc_count = 0
    
    for rec in records:
        seq = rec.sequence.upper()
        contigs += 1
        total_length += len(seq)
        gc_count += seq.count('G') + seq.count('C')
    
    gc_content = (100.0 * gc_count / total_length) if total_length > 0 else None
    
    logging.info(f"Loaded FASTA with {contigs} contigs, {total_length:,} bp")
    
    return {
        'denovo_contigs': contigs,
        'denovo_length': total_length,
        'denovo_gc': gc_content
    }


def calculate_taxid_stats(
    data: Dict[str, Any],
    ncbi_stats: List[Dict[str, str]],
    taxid: Optional[str],
    min_n: int
) -> Dict[str, Any]:
    """
    Calculate z-scores for assembly metrics using NCBI reference stats.
    
    Matches by taxid if provided, otherwise by species name.
    Only calculates z-scores if reference has sufficient samples (>= min_n).
    """
    result: Dict[str, Any] = {}

    # Determine what to match on
    target = taxid if taxid else data.get('species')
    match_column = 'taxids' if taxid else 'names'

    target = normalize_value(target)
    if not target:
        return result
            
    # Find matching reference stats
    for rec in ncbi_stats:
        rec_values = [ normalize_value(v) for v in rec.get(match_column, []) ]
        
        if target not in rec_values:
            continue
        
        # Extract reference statistics
        length_mean = float(rec.get("length_mean", 0))
        length_sd = float(rec["length_stdev"]) if rec.get("length_stdev") is not None else None
        gc_mean = float(rec.get("gc_mean", 0))
        gc_sd = float(rec["gc_stdev"]) if rec.get("gc_stdev") is not None else None
        n_samples = int(rec.get("n", 0))

        # Estimate sequencing depth
        total_bases = data.get('total_bases_after_filtering', 0)
        if length_mean > 0:
            result['estimated_depth'] = int(round(float(total_bases) / float(length_mean)))
        else:
            result['estimated_depth'] = 0

        # Only calculate z-scores if sufficient reference samples
        if n_samples < min_n:
            logging.warning(f"Insufficient reference samples to calculate z-scores (n={n_samples}, min={min_n})")
            return result
    
        # Calculate z-scores
        if 'denovo_length' not in data:
            continue

        assembly_length = float(data.get('denovo_length', 0))
        assembly_gc = float(data.get('denovo_gc', 0))
        
        result['denovo_length_z'] = calculate_z_score(
            assembly_length, length_mean, length_sd
        )
        result['denovo_gc_z'] = calculate_z_score(
            assembly_gc, gc_mean, gc_sd
        )

        logging.info(f"Calculated z-scores using {n_samples} reference genomes")
        break

    return result


def parse_read_stats(stats_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Extract read statistics from fastp summary."""
    result: Dict[str, Any] = {}
    
    summary = stats_dict.get('summary', {})
    for stage in ['before_filtering', 'after_filtering']:
        if stage in summary:
            for key, value in summary[stage].items():
                result[f"{key}_{stage}"] = value
    
    return result


def parse_species(species_dict: Dict[str, str]) -> tuple[Dict[str, Any], Optional[str]]:
    """Extract species identification from GAMBIT output."""
    species_name = species_dict.get('predicted.name')
    rank = species_dict.get('predicted.rank')
    taxid = species_dict.get('predicted.ncbi_id')

    # Fall back to next best match if prediction unavailable
    if not species_name:
        species_name = species_dict.get('next.name')

    result = {'species': species_name}
    return result, taxid


def parse_subtype(subtype_dict: Dict[str, str]) -> Dict[str, Any]:
    """Extract subtype information from subtyping results."""
    subtype_name = subtype_dict.get('subtype')
    confidence = subtype_dict.get('closest_ani')
    
    # Extract base subtype (remove suffix after hyphen)
    if subtype_name:
        subtype_name = subtype_name.split('-', 1)[0]
    
    return {
        'subtype': subtype_name,
        'subtype_ani': confidence
    }


def parse_samplesheet(
    samplesheet_data: List[Dict[str, str]],
    sample_name: str
) -> Dict[str, Any]:
    """Extract sample-specific data from samplesheet."""
    for record in samplesheet_data:
        if record.get('sample') == sample_name:
            logging.info(f"Found sample in samplesheet: {list(record.keys())}")
            return record
    
    logging.warning(f"Sample '{sample_name}' not found in samplesheet")
    return {}


# ----------------------------
# QC Functions
# ----------------------------

def perform_auto_qc(
    data: Dict[str, Any],
    min_depth: int,
    min_qual: float,
    max_z: float
) -> Dict[str, Any]:
    """
    Perform automated quality control checks.
    
    Checks:
    - Species identified
    - Subtype identified
    - Q30 rate >= threshold
    - Read depth >= threshold
    - Assembly length z-score < threshold (if available)
    - GC content z-score < threshold (if available)
    """
    qc_status = 'PASS'
    qc_und: List[str] = []
    qc_fail: List[str] = []
    qc_error: List[str] = []


    # Define QC criteria (field: (operator, threshold) or None for required)
    qc_criteria = {
        'q30_rate_after_filtering': ('>=', min_qual),
        'species': None,
        'subtype': None,
        'estimated_depth': ('>=', min_depth),
        'denovo_length_z': ('<', max_z),
        'denovo_gc_z': ('<', max_z)
    }
    
    for field, criterion in qc_criteria.items():
        value = data.get(field)

        # Check if field exists
        if value is None:
            qc_und.append(field)
            if field.startswith('denovo'):
                continue
            qc_status = 'FAIL'
            break

        # If no criterion specified, just check for presence
        if criterion is None:
            continue

        # Compare against threshold
        try:
            operator, threshold = criterion
            label = field

            # For z-scores, check absolute value
            if field.endswith('_z'):
                value = abs(value)
                label = f"|{field}|"

            if not compare_values(value, operator, threshold):
                qc_status = 'FAIL'
                observed = round(value, 2) if isinstance(value, float) else value
                qc_fail.append(f"{label} = {observed} (required {operator} {threshold})")
        except Exception as e:
            qc_status = 'FAIL'
            qc_error.append(field)
            logging.error(f"QC comparison failed for {field}: {e}")
    
    qc_reasons = [ 
        f"Undetermined: {', '.join(qc_und)}" if qc_und else '', 
        f"Failure: {', '.join(qc_fail)}" if qc_fail else '', 
        f"Error: {', '.join(qc_error)}" if qc_error else ''
    ]

    data['qc_status'] = qc_status
    data['qc_reason'] = '; '.join([r for r in qc_reasons if r])
    
    logging.info(f"QC Status: {data['qc_status']}")
    if qc_reasons:
        logging.info(f"QC Reasons: {data['qc_reason']}")
    
    return data


# ----------------------------
# Main Function
# ----------------------------

def main():
    """Main workflow summarization function."""
    VERSION = "1.2"

    parser = argparse.ArgumentParser(
        description="Summarize outputs from bioinformatics workflows",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--sample", required=True,
                        help="Sample name")
    parser.add_argument("--samplesheet", required=True,
                        help="Samplesheet CSV path")
    
    # Optional input files
    parser.add_argument("--ncbi_stats",
                        help="NCBI reference statistics CSV")
    parser.add_argument("--read_stats",
                        help="Fastp summary output (JSON)")
    parser.add_argument("--species",
                        help="GAMBIT species identification output (CSV)")
    parser.add_argument("--subtype",
                        help="Subtype identification output (CSV)")
    parser.add_argument("--denovo",
                        help="De novo assembly (FASTA)")
    
    # QC parameters
    parser.add_argument("--min_ncbi_stats_n", type=int, default=3,
                        help="Minimum samples in NCBI stats for z-score calculation")
    parser.add_argument("--min_depth", type=int, default=30,
                        help="Minimum read depth for QC pass")
    parser.add_argument("--min_qual", type=float, default=0.8,
                        help="Minimum Q30 rate for QC pass")
    parser.add_argument("--max_z_score", type=float, default=2.58,
                        help="Maximum absolute z-score for QC pass")
    
    # Logging options
    parser.add_argument("--log-level",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO',
                        help="Logging verbosity level")
    parser.add_argument("--log-file",
                        help="Optional log file path")
    
    parser.add_argument("--version", action="version",
                        version=VERSION)

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        filename=args.log_file,
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    logging.info(f"Starting workflow summary for sample: {args.sample}")

    # Initialize result dictionary
    data: Dict[str, Any] = {'sample': args.sample}

    # Parse samplesheet (for manual overrides)
    samplesheet_data = load_csv(args.samplesheet, 'samplesheet')
    data.update(parse_samplesheet(samplesheet_data, args.sample))

    # Parse species identification
    taxid: Optional[str] = None
    if data.get('species'):
        logging.info("Using user-supplied species from samplesheet")
    elif args.species:
        logging.info("Extracting species identification")
        species_data = load_csv(args.species, 'species')
        if species_data:
            species_parsed, taxid = parse_species(species_data[0])
            data.update(species_parsed)
        else:
            logging.warning("Species CSV is empty")
    else:
        logging.warning("No species data provided")

    # Parse subtype
    if data.get('subtype'):
        logging.info("Using user-supplied subtype from samplesheet")
    elif args.subtype:
        logging.info("Extracting subtype information")
        subtype_data = load_csv(args.subtype, 'subtype')
        if subtype_data:
            subtype_parsed = parse_subtype(subtype_data[0])
            data.update(subtype_parsed)
        else:
            logging.warning("Subtype CSV is empty")
    else:
        logging.warning("No subtype data provided")

    # Parse assembly statistics
    if args.denovo:
        logging.info("Calculating genome statistics from assembly")
        with load_fasta(args.denovo) as fasta_records:
            genome_stats_data = calculate_genome_stats(fasta_records)
        data.update(genome_stats_data)
    else:
        logging.warning("No assembly data provided")

    # Parse read statistics
    if args.read_stats:
        logging.info("Extracting read statistics")
        stats = load_json(args.read_stats, 'read_stats')
        stats_parsed = parse_read_stats(stats)
        data.update(stats_parsed)
    else:
        logging.warning("No read statistics provided")

    # Calculate z-scores using NCBI reference data
    if args.ncbi_stats:
        logging.info("Calculating z-scores from NCBI reference statistics")
        ncbi_stats_data = load_json(args.ncbi_stats, 'ncbi_stats')
        taxid_stats_data = calculate_taxid_stats(
            data, ncbi_stats_data, taxid, args.min_ncbi_stats_n
        )
        if taxid_stats_data:
            data.update(taxid_stats_data)
        else:
            logging.warning("Could not calculate NCBI-based statistics")
    else:
        logging.warning("NCBI reference statistics not provided")

    # Perform automated QC
    data = perform_auto_qc(data, args.min_depth, args.min_qual, args.max_z_score)

    # Format output values
    for key, value in list(data.items()):
        # Convert lists to semicolon-delimited strings
        if isinstance(value, list):
            data[key] = ';'.join(map(str, value))
        # Round floats (except depth fields which should be integers)
        elif isinstance(value, float) and not key.endswith('_depth'):
            data[key] = round(value, 2)

    # Samples processed here always come from the input samplesheet
    data['status'] = 'new'

    # Define output column order
    output_columns = [
        'sample',
        'status',
        'qc_status',
        'qc_reason',
        'species',
        'subtype',
        'subtype_ani',
        'estimated_depth',
        'denovo_contigs',
        'denovo_length',
        'denovo_length_z',
        'denovo_gc',
        'denovo_gc_z',
        'total_reads_after_filtering',
        'total_bases_after_filtering',
        'q30_bases_after_filtering',
        'q30_rate_after_filtering',
        'read1_mean_length_after_filtering',
        'read2_mean_length_after_filtering',
        'gc_content_after_filtering',
        'total_reads_before_filtering',
        'total_bases_before_filtering',
        'q30_bases_before_filtering',
        'q30_rate_before_filtering',
        'read1_mean_length_before_filtering',
        'read2_mean_length_before_filtering',
        'gc_content_before_filtering',
    ]
    
    # Filter data to only include specified columns (in order)
    filtered_data = {col: data.get(col, '') for col in output_columns}

    # Write output CSV
    output_file = f"{args.sample}_summary.csv"
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=output_columns)
        writer.writeheader()
        writer.writerow(filtered_data)

    logging.info(f"Summary written to: {output_file}")


if __name__ == '__main__':
    main()