#!/usr/bin/env python3
"""
Summarize outputs from bioinformatics workflows.

Aggregates read statistics, species identification, assembly metrics,
and performs automated QC checks.
"""

import json
import csv
import re
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


def sanitize(value: Any) -> Optional[str]:
    """Normalize a name the way the pipeline does (Utils.sanitize)."""
    if value is None:
        return None
    return re.sub(r'[^A-Za-z0-9_-]', '_', str(value).strip()).lower()


def parse_range(value: Any) -> Optional[tuple]:
    """Return (min, max) floats from a [min, max] list, or None if invalid."""
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        return None
    try:
        low, high = float(value[0]), float(value[1])
    except (TypeError, ValueError):
        return None
    return (low, high) if low <= high else None


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


def find_ncbi_record(
    ncbi_stats: List[Dict[str, Any]],
    taxid: Optional[str],
    species: Optional[str],
    aliases: Optional[List[str]] = None
) -> Optional[Dict[str, Any]]:
    """
    Find the NCBI stats record for a sample: by taxid if provided, then by
    species name, then by the species' aliases from the reference manifest
    (e.g. 'Candida auris' -> 'Candidozyma auris').
    """
    targets = [(taxid, 'taxids'), (species, 'names')] + [(a, 'names') for a in (aliases or [])]
    for target, column in targets:
        target = normalize_value(target)
        if not target:
            continue
        for rec in ncbi_stats:
            if target in [normalize_value(v) for v in rec.get(column, [])]:
                return rec
    return None


def find_reference_record(
    reference_qc: List[Dict[str, Any]],
    species: Optional[str]
) -> Optional[Dict[str, Any]]:
    """Find the reference manifest QC entry whose species names include the sample's species."""
    target = sanitize(species)
    if not target:
        return None
    for rec in reference_qc:
        if target in [sanitize(v) for v in rec.get('species', [])]:
            return rec
    return None


def resolve_qc_ranges(
    ncbi_rec: Optional[Dict[str, Any]],
    ref_rec: Optional[Dict[str, Any]],
    min_n: int
) -> Dict[str, Dict[str, Any]]:
    """
    Acceptable assembly length and GC ranges, per metric: from the reference
    manifest when it sets one, otherwise from NCBI (if it has >= min_n genomes).

    Returns e.g. {'length': {'range': (min, max), 'source': 'manifest'}, ...}
    """
    ranges: Dict[str, Dict[str, Any]] = {}
    for metric in ('length', 'gc'):
        key = f'{metric}_range'

        rng = parse_range(ref_rec.get(key)) if ref_rec else None
        if rng:
            ranges[metric] = {'range': rng, 'source': 'manifest'}
            continue

        if not ncbi_rec:
            continue
        n_samples = int(ncbi_rec.get('n', 0))
        if n_samples < min_n:
            logging.warning(f"Insufficient NCBI genomes for the {metric} range (n={n_samples}, min={min_n})")
            continue
        rng = parse_range(ncbi_rec.get(key))
        if rng:
            ranges[metric] = {'range': rng, 'source': 'ncbi'}
        else:
            logging.warning(f"NCBI stats have no '{key}'; regenerate them with make_ncbi_stats.py")

    for metric, r in ranges.items():
        logging.info(f"{metric} range {r['range']} from {r['source']}")
    return ranges


def estimate_depth(
    data: Dict[str, Any],
    ncbi_rec: Optional[Dict[str, Any]],
    ranges: Dict[str, Dict[str, Any]]
) -> Optional[int]:
    """
    Estimated sequencing depth: filtered bases / expected genome length. The
    expected length is the NCBI mean length, or the midpoint of the manifest
    length range for species without NCBI data. None if neither is available.
    """
    genome_length = float(ncbi_rec.get('length_mean') or 0) if ncbi_rec else 0
    if genome_length <= 0 and ranges.get('length', {}).get('source') == 'manifest':
        low, high = ranges['length']['range']
        genome_length = (low + high) / 2
    if genome_length <= 0:
        return None
    total_bases = float(data.get('total_bases_after_filtering') or 0)
    return int(round(total_bases / genome_length))


def format_range(rng: tuple, digits: int) -> str:
    """Format (min, max) as 'min-max'."""
    return f"{rng[0]:.{digits}f}-{rng[1]:.{digits}f}"


def describe_sources(ranges: Dict[str, Dict[str, Any]]) -> str:
    """e.g. 'manifest', or 'length: manifest; gc: ncbi' when they differ."""
    sources = {m: r['source'] for m, r in ranges.items()}
    if not sources:
        return ''
    if len(set(sources.values())) == 1 and len(sources) == 2:
        return next(iter(sources.values()))
    return '; '.join(f"{m}: {src}" for m, src in sources.items())


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
    ranges: Dict[str, Dict[str, Any]],
    classify_fail: Optional[str] = None
) -> Dict[str, Any]:
    """
    Perform automated quality control checks.

    A sample the pipeline could not match to a reference (classify_fail set)
    always fails, with that reason reported first.

    Checks:
    - Species identified
    - Subtype identified
    - Q30 rate >= threshold
    - Read depth >= threshold
    - Assembly length within the species' range (if available)
    - Assembly GC content within the species' range (if available)
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
    }

    missing_required = False
    for field, criterion in qc_criteria.items():
        value = data.get(field)

        # Missing required fields fail QC and stop further checks
        if value is None:
            qc_und.append(field)
            qc_status = 'FAIL'
            missing_required = True
            break

        # If no criterion specified, just check for presence
        if criterion is None:
            continue

        # Compare against threshold
        try:
            operator, threshold = criterion
            if not compare_values(value, operator, threshold):
                qc_status = 'FAIL'
                observed = round(value, 2) if isinstance(value, float) else value
                qc_fail.append(f"{field} = {observed} (required {operator} {threshold})")
        except Exception as e:
            qc_status = 'FAIL'
            qc_error.append(field)
            logging.error(f"QC comparison failed for {field}: {e}")

    # Assembly length / GC must fall within the species' range. A missing
    # assembly or range is reported as undetermined but does not fail QC.
    if not missing_required:
        for metric, field, digits in (('length', 'denovo_length', 0), ('gc', 'denovo_gc', 2)):
            value = data.get(field)
            rng = ranges.get(metric, {}).get('range')
            if value is None or rng is None:
                qc_und.append(f"{field}_range" if value is not None else field)
                continue
            try:
                low, high = rng
                if not low <= float(value) <= high:
                    qc_status = 'FAIL'
                    observed = round(value, 2) if isinstance(value, float) else value
                    qc_fail.append(f"{field} = {observed} (required {format_range(rng, digits)})")
            except Exception as e:
                qc_status = 'FAIL'
                qc_error.append(field)
                logging.error(f"QC comparison failed for {field}: {e}")

    if classify_fail:
        qc_status = 'FAIL'

    qc_reasons = [
        f"Classification: {classify_fail}" if classify_fail else '',
        f"Undetermined: {', '.join(qc_und)}" if qc_und else '',
        f"Failure: {', '.join(qc_fail)}" if qc_fail else '',
        f"Error: {', '.join(qc_error)}" if qc_error else ''
    ]

    data['qc_status'] = qc_status
    data['qc_reason'] = '; '.join([r for r in qc_reasons if r])

    logging.info(f"QC Status: {data['qc_status']}")
    if data['qc_reason']:
        logging.info(f"QC Reasons: {data['qc_reason']}")

    return data


# ----------------------------
# Main Function
# ----------------------------

def main():
    """Main workflow summarization function."""
    VERSION = "1.4"

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
                        help="NCBI species statistics (JSON, from make_ncbi_stats.py)")
    parser.add_argument("--reference_qc",
                        help="Species QC ranges from the reference manifest (JSON); "
                             "used instead of the NCBI ranges where set")
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
                        help="Minimum NCBI genomes for the NCBI length / GC ranges to be used")
    parser.add_argument("--min_depth", type=int, default=30,
                        help="Minimum read depth for QC pass")
    parser.add_argument("--min_qual", type=float, default=0.8,
                        help="Minimum Q30 rate for QC pass")
    parser.add_argument("--classify_fail",
                        help="Reason the sample could not be matched to a reference; "
                             "forces qc_status FAIL")
    
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

    # Expected assembly length / GC ranges (manifest over NCBI) and depth
    ref_rec = None
    if args.reference_qc:
        ref_rec = find_reference_record(load_json(args.reference_qc, 'reference_qc'), data.get('species'))

    ncbi_rec = None
    if args.ncbi_stats:
        aliases = ref_rec.get('species', []) if ref_rec else []
        ncbi_rec = find_ncbi_record(load_json(args.ncbi_stats, 'ncbi_stats'), taxid, data.get('species'), aliases)
        if not ncbi_rec:
            logging.warning("Species not found in NCBI statistics")
    else:
        logging.warning("NCBI statistics not provided")

    ranges = resolve_qc_ranges(ncbi_rec, ref_rec, args.min_ncbi_stats_n)
    data['estimated_depth'] = estimate_depth(data, ncbi_rec, ranges)
    if 'length' in ranges:
        data['denovo_length_range'] = format_range(ranges['length']['range'], 0)
    if 'gc' in ranges:
        data['denovo_gc_range'] = format_range(ranges['gc']['range'], 2)
    data['qc_range_source'] = describe_sources(ranges)

    # Perform automated QC
    data = perform_auto_qc(data, args.min_depth, args.min_qual, ranges, args.classify_fail)

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
        'denovo_length_range',
        'denovo_gc',
        'denovo_gc_range',
        'qc_range_source',
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