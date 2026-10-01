#!/usr/bin/env python3

import json
import csv
import screed
import argparse
import logging
from typing import List, Dict, Any, Generator


# ----------------------------
# File helper functions
# ----------------------------

def load_csv(path: str, source: str = None) -> List[Dict[str, str]]:
    """Load CSV into a list of dictionaries (rows)."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        data = [row for row in reader]
    logging.info(f"Loaded CSV with {len(data)} rows ({source})")
    return data
    
def parse_amr(amr_dict):
    res = {'amr_variants_target': [], 'amr_variants_other': []}
    for rec in amr_dict:        
        gene   = rec.get('Target_Gene', None)
        region = rec.get('Target_Region', None)
        impact = rec.get('Annotation_Impact', None)
        if not gene or impact not in ['MODERATE', 'HIGH']:
            continue
        mutation = rec.get('HGVS.p')
        if region:
            res['amr_variants_target'].append(f"{gene}({region}):{mutation}")
        else:
            res['amr_variants_other'].append(f"{gene}:{mutation}")

    return { k: ';'.join(v) for k, v in res.items()}

# ----------------------------
# Main
# ----------------------------

def main():
    version = "1.1"

    parser = argparse.ArgumentParser(
        description="Summarize outputs from various workflows"
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--summaryline", required=True, help="Summaryline path")
    parser.add_argument("--amr", required=True, help="snpeff_parser.py summary output (CSV)")
    parser.add_argument("--log-level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help="Logging level")
    parser.add_argument("--version", action="version", version=version, help="Show script version and exit.")
    parser.add_argument("--log-file", help="Log file path (optional)")
    
    args = parser.parse_args()

    # configure logging
    logging.basicConfig(
        filename=args.log_file,
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    summaryline_data = load_csv(args.summaryline, 'summaryline')
    summaryline_data = summaryline_data[0]

    amr_data = load_csv(args.amr, 'amr')
    amr_parsed = parse_amr(amr_data)

    # Insert AMR columns after the assembly QC columns ('denovo_gc_z' in
    # summaries written before QC ranges replaced z-scores)
    anchor = next((k for k in ('qc_range_source', 'denovo_gc_z') if k in summaryline_data), None)
    final = {}
    if anchor:
        for k, v in summaryline_data.items():
            final[k] = v
            if k == anchor:
                final = final | amr_parsed
    else:
        final = summaryline_data | amr_parsed

    with open(f"{args.sample}_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(final.keys()))
        writer.writeheader()
        writer.writerow(final)


if __name__ == '__main__':
    main()
