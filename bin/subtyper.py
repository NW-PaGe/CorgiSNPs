#!/usr/bin/env python3

"""
subtyper.py

Subtype assemblies against reference assemblies using sourmash.

Two commands:

    sketch  Build one signature per reference, with the subtype stored in the
            signature name. The subtype TSV maps reference paths to subtypes.

                subtyper.py sketch subtypes.tsv ref1.fa ref2.fa ref3.fa \\
                    --prefix candidozyma_auris

            subtypes.tsv columns: subtype, path

    assign  Compare an assembly to those signatures. The subtype is the name
            of the signature with the highest ANI to the sample, provided that
            ANI meets the run-wide threshold; otherwise it is 'undefined'.

                subtyper.py assign candidozyma_auris.sig.zip contigs.fa \\
                    --threshold 0.997

Authors:
    Jared Johnson, jared.johnson@doh.wa.gov
    Zack Mudge, ZMudge@cdc.gov
"""

import argparse
import csv
import logging
import os
import sys
from typing import Any, Dict, List, Tuple

import screed
import sourmash as sm
from sourmash.save_load import SaveSignaturesToLocation

VERSION = "v3.1.0"


# ----------------------------

# ----------------------------
def setup_logging(log_level: str = 'INFO', log_file: str | None = None) -> logging.Logger:
    level = getattr(logging, log_level.upper())
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(level) 
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(level)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


# ----------------------------
# Shared helpers
# ----------------------------
def parse_threshold(raw: float) -> float:
    """Convert an ANI threshold to a percentage (0.997 -> 99.7; 99.7 stays 99.7)."""
    try:
        value = float(raw)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid ANI threshold: {raw}") from e

    if raw > 1.0 or raw < 0.0:
        raise ValueError(f"ANI threshold must be 0-1: {raw}")

    return raw * 100.0


def read_tsv(path: str, required: set) -> List[Dict[str, str]]:
    """Read a TSV, checking that the required columns are present."""
    with open(path, newline='', encoding='utf-8') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        missing = required - set(rdr.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(sorted(missing))}")
        return [ {k: (v or '').strip() for k, v in row.items()} for row in rdr ]


def sketch_seq(filepath: str, ksize: int, scaled: int) -> Tuple[sm.MinHash, int, int]:
    """Build a MinHash from a sequence file (FASTA/FASTQ, optionally gzipped), adding each record separately."""
    mh = sm.MinHash(n=0, ksize=ksize, scaled=scaled)
    records = 0
    bases = 0
    for rec in screed.open(filepath):
        mh.add_sequence(rec.sequence, force=True)
        records += 1
        bases += len(rec.sequence)
    if not records:
        raise ValueError(f"No sequences found in {filepath}")
    logging.info(f"Sketched {os.path.basename(filepath)}: {records:,} records; {bases:,} bp; {len(mh):,} hashes")
    return mh, records, bases


def reference_subtypes(path: str) -> Dict[str, str]:
    """Read reference-to-subtype mappings and index them by path and basename."""
    rows = [r for r in read_tsv(path, {'subtype', 'path'}) if r['subtype'] and r['path']]
    if not rows:
        raise ValueError(f"No reference mappings listed in {path}")

    mappings: Dict[str, str] = {}
    for row in rows:
        keys = {os.path.normpath(row['path']), os.path.basename(row['path'])}
        for key in keys:
            previous = mappings.get(key)
            if previous is not None and previous != row['subtype']:
                raise ValueError(
                    f"Reference '{key}' maps to multiple subtypes in {path}: "
                    f"'{previous}' and '{row['subtype']}'"
                )
            mappings[key] = row['subtype']
    return mappings


def sample_name(path: str) -> str:
    """Return a sample name from a FASTA/FASTQ path, removing common suffixes."""
    name = os.path.basename(path)
    if name.endswith('.gz'):
        name = name[:-3]
    for suffix in ('.fasta', '.fastq', '.fa', '.fq', '.fna'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]


# ----------------------------
# sketch
# ----------------------------
def run_sketch(args: argparse.Namespace) -> None:
    """Create one signature per supplied reference, with subtype stored in its name."""
    mappings = reference_subtypes(args.subtypes)
    references = []
    for path in args.references:
        if not os.path.isfile(path):
            raise ValueError(f"Reference assembly does not exist: {path}")

        subtype = mappings.get(os.path.normpath(path)) or mappings.get(os.path.basename(path))
        if subtype is None:
            raise ValueError(f"No subtype mapping in {args.subtypes} for reference: {path}")
        references.append((path, subtype))

    sig_out = f"{args.prefix}.sig.zip"
    logging.info(f"Sketching {len(references)} reference(s); ksize={args.ksize}; scaled={args.scaled}")

    with SaveSignaturesToLocation(sig_out) as save_sigs:
        for path, subtype in references:
            mh, _, _ = sketch_seq(path, args.ksize, args.scaled)
            save_sigs.add(sm.SourmashSignature(mh, name=subtype, filename=os.path.basename(path)))

    logging.info(f"Wrote: {sig_out}")


# ----------------------------
# assign
# ----------------------------
def load_signatures(sig_path: str) -> Tuple[List[Any], int, int]:
    """Load all signatures; ensure a single ksize and scaled, and that every signature is named."""
    sigs = list(sm.load_file_as_signatures(sig_path))
    if not sigs:
        raise ValueError(f"No signatures found in {sig_path}")

    kset = {s.minhash.ksize for s in sigs}
    sset = {s.minhash.scaled for s in sigs}
    if len(kset) > 1:
        raise ValueError(f"Multiple ksize values in {sig_path}: {kset}")
    if len(sset) > 1:
        raise ValueError(f"Multiple scaled values in {sig_path}: {sset}")
    if any(not s.name for s in sigs):
        raise ValueError(f"Signature(s) in {sig_path} have no name; names are used as subtypes")

    ksize, scaled = next(iter(kset)), next(iter(sset))
    logging.info(f"Reference signatures: {len(sigs)}; ksize={ksize}; scaled={scaled}")
    return sigs, ksize, scaled


def write_assign_csv(path: str, row: dict) -> None:
    """Write a single-row CSV (overwrite)."""
    fields = [
        'sample', 'signature_file', 'ksize', 'scaled', 'n_db_signatures',
        'subtype', 'closest_subtype', 'closest_ani', 'threshold', 'passed_threshold'
    ]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerow(row)


def run_assign(args: argparse.Namespace) -> None:
    """Assign a subtype to one assembly using reference signatures."""
    sample = args.sample or sample_name(args.assembly)
    out_csv = args.out or f"{sample}_subtype.csv"

    sigs, ksize, scaled = load_signatures(args.signatures)
    threshold = parse_threshold(args.threshold)
    sample_mh, _, _ = sketch_seq(args.assembly, ksize, scaled)

    # Compare each reference signature to the sample sketch
    top_ani = 0.0
    top_matches: set[str] = set()

    for sig in sigs:
        res = sig.minhash.containment_ani(sample_mh)
        try:
            dist = res.dist
        except AttributeError:
            dist = float(res)
        ani = 100.0 * (1.0 - dist) if dist is not None else 0.0
        logging.debug(f"{sig.name}: ANI {ani:.4f}")

        if ani > top_ani:
            top_ani = ani
            top_matches = { sig.name }
        elif ani == top_ani and ani > 0:
            top_matches.add(sig.name)

    closest_subtype = " / ".join(sorted(top_matches))

    passed = bool(top_matches) and top_ani >= threshold
    subtype = closest_subtype if passed else 'undefined'

    if top_matches:
        logging.info(f"Closest subtype: {closest_subtype} (ANI {top_ani:.2f}%, threshold {threshold}%)")
    else:
        logging.warning("No reference signature shared any k-mers with the sample")
    logging.info(f"Assigned subtype: {subtype}")

    write_assign_csv(out_csv, {
        'sample': sample,
        'signature_file': args.signatures,
        'ksize': ksize,
        'scaled': scaled,
        'n_db_signatures': len(sigs),
        'subtype': subtype,
        'closest_subtype': closest_subtype,
        'closest_ani': round(top_ani, 2) if top_matches else '',
        'threshold': threshold,
        'passed_threshold': passed
    })
    logging.info(f"Wrote: {out_csv}")


# ----------------------------
# Main
# ----------------------------
def main():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--log-level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help="Logging level")
    common.add_argument("--log-file", help="Log file path (optional)")

    p = argparse.ArgumentParser(description='Subtype assemblies against reference signatures with sourmash')
    p.add_argument("--version", action="version", version=VERSION)
    sub = p.add_subparsers(dest='command', required=True)

    ps = sub.add_parser('sketch', parents=[common],
                        help='Create reference signatures named by subtype')
    ps.add_argument('subtypes', help='TSV mapping reference paths to subtypes (columns: subtype, path)')
    ps.add_argument('references', nargs='+', help='Reference assembly path(s)')
    ps.add_argument('--prefix', required=True,
                    help='Output prefix: writes <prefix>.sig.zip')
    ps.add_argument('--ksize', type=int, default=31, help='k-mer size (default: 31)')
    ps.add_argument('--scaled', type=int, default=100, help='Scaled value (default: 100)')
    ps.set_defaults(func=run_sketch)

    pa = sub.add_parser('assign', parents=[common],
                        help='Assign a subtype to an assembly using reference signatures')
    pa.add_argument('signatures', help='Signatures from the sketch command')
    pa.add_argument('assembly', help='Assembly or reads to subtype (FASTA/FASTQ)')
    pa.add_argument('--threshold', type=float, default=0.997,
                    help='ANI threshold as a fraction (default: 0.997)')
    pa.add_argument('--sample', help='Sample name (default: derived from assembly filename)')
    pa.add_argument('--out', help='Output CSV path (default: <sample>_subtype.csv)')
    pa.set_defaults(func=run_assign)

    args = p.parse_args()
    setup_logging(args.log_level, args.log_file)
    args.func(args)


if __name__ == '__main__':
    main()
