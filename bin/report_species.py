#!/usr/bin/env python3

import re
import json
import csv
import gzip
import argparse
import logging
from typing import List, Dict, Any, TextIO
from Bio import Phylo
import numpy as np
from sklearn.cluster import DBSCAN
from pathlib import Path
from datetime import datetime, timezone

# ----------------------------
# Microreact helpers
# ----------------------------

def _ensure_path(p: Path, label: str) -> bool:
    if not p.exists():
        logging.warning(f"{label} missing: {p}")
        return False
    return True

def _attach_text_file(mr: Dict[str, Any], slot: str, path: Path, outname: str) -> None:
    mr.setdefault('files', {}).setdefault(slot, {})
    mr['files'][slot]['blob'] = path.read_text(encoding='utf-8')
    mr['files'][slot]['name'] = outname

# ----------------------------
# Naming helpers
# ----------------------------

def make_prefix(species: str, subtype: str) -> str:
    """
    Build the file prefix from species and subtype: joined with a dash,
    lowercased, and every whitespace character replaced by an underscore.
    e.g. ("Escherichia coli", "ST 131") -> "escherichia_coli-st_131"
    """
    joined = f"{species.strip()}-{subtype.strip()}".lower()
    return re.sub(r'\s', '_', joined)

def make_meta_name(species: str, subtype: str, date: str) -> str:
    """Microreact project name: '<species> <subtype> (<date>)'."""
    return f"{species.strip()} {subtype.strip()} ({date})"

# ----------------------------
# File helpers
# ----------------------------

def open_text(path: str) -> TextIO:
    """
    Open a file for reading as text. Gzipped files are decompressed on the
    fly; they are detected by their content, so the '.gz' extension is not
    required.
    """
    with open(path, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'
    if is_gzip:
        return gzip.open(path, 'rt', encoding='utf-8', newline='')
    return open(path, 'r', encoding='utf-8', newline='')

def detect_sep(path: str) -> str:
    """Tab for .tsv / .tsv.gz files, comma otherwise."""
    name = path.lower()
    if name.endswith('.gz'):
        name = name[:-3]
    return '\t' if name.endswith('.tsv') else ','

def load_json(path: str) -> Dict[str, Any]:
    with open_text(path) as f:
        return json.load(f)

def load_csv(path: str, sep: str = ',') -> List[Dict[str, str]]:
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter=sep))

def load_dist(path: str, out_matrix: str = 'matrix.csv'):
    rowids, colids, mr_out = [], None, []
    M = None
    with open_text(path) as f:
        for rn, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            if rn == 0:
                colids = line.split(',')[1:]
                n = len(colids)
                M = np.empty((n, n), dtype=float)
                mr_out.append(['sample_id']+colids)
                continue
            row = line.split(',')
            id_row, values = row[0], row[1:]
            rowids.append(id_row)
            for cn, _ in enumerate(colids):
                M[rn - 1, cn] = float(values[cn])
            mr_out.append(row)
    if rowids != colids:
        raise ValueError(f"Row/column IDs mismatch in {path}")
    with open(out_matrix, 'w', encoding='utf-8') as f:
        f.write('\n'.join([','.join(r) for r in mr_out]))
    return colids, M

def list2map(data: List[Dict[str, Any]], key: str) -> Dict[str, Dict[str, Any]]:
    return {row[key]: row for row in data if key in row}

def patristic_distance_matrix(tree):
    tips = tree.get_terminals()
    labels = [t.name for t in tips]
    n = len(tips)
    M = np.zeros((n, n), dtype=float)
    for i, a in enumerate(tips):
        for j in range(i, n):
            d = tree.distance(a, tips[j])
            M[i, j] = M[j, i] = d
    return labels, M

def partition_tree(D, eps, min_samples=1):
    db = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    y = db.fit_predict(D)
    max_cluster = max(y) if len(y) else -1
    clusters = [c + 1 for c in y]
    for i, c in enumerate(clusters):
        if c == 0:
            c = max_cluster + 2
            max_cluster = c - 1
        clusters[i] = c
    return clusters

def find_links(D, labels, strong_thresh, inter_thresh):
    res = {l: {'strong_links': [], 'inter_links': []} for l in labels}
    n = len(labels)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = D[i, j]
            if d <= strong_thresh:
                res[labels[i]]['strong_links'].append(labels[j])
            elif d <= inter_thresh:
                res[labels[i]]['inter_links'].append(labels[j])
    return res

def extend_dict(main: Dict[str, Dict[str, Any]], new: Dict[str, Dict[str, Any]]):
    for k, v in new.items():
        main[k] = main.get(k, {}) | v

def is_reference_sample(sample_id: str) -> bool:
    return bool(sample_id) and (sample_id == 'Reference' or sample_id.startswith('Reference_'))

def reference_row(stats: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    if 'Reference' in stats:
        return stats['Reference']
    for key, row in stats.items():
        if key.startswith('Reference_'):
            return row
    return {}

# ----------------------------
# Main
# ----------------------------

def main():
    VERSION = "1.4"

    parser = argparse.ArgumentParser(description="Summarize outputs from various workflows. All inputs may be gzipped.")
    parser.add_argument("--species", required=True, help="Species name (e.g. 'Escherichia coli').")
    parser.add_argument("--subtype", required=True, help="Subtype (e.g. 'ST131').")
    parser.add_argument("--aln_stats", required=True, help="Core alignment stats from PolyCore.")
    parser.add_argument("--summary", required=True, help="Combined summary file. May contain more than what is in core alignment.")
    parser.add_argument("--tree")
    parser.add_argument("--dist")
    parser.add_argument("--microreact")
    parser.add_argument("--partition_distance", default=100, type=float)
    parser.add_argument("--strong_link", default=5, type=float)
    parser.add_argument("--inter_link", default=10, type=float)
    parser.add_argument("--log-level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO')
    parser.add_argument("--log-file")
    parser.add_argument("--version", action="version",
                        version=VERSION)
    args = parser.parse_args()

    logging.basicConfig(
        filename=args.log_file,
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    logging.info("Starting")

    # Single run time, used for both the meta name date and the timestamp
    run_time = datetime.now(timezone.utc)

    # Accumulators
    data: Dict[str, Dict[str, Any]] = {}
    stats_dict: Dict[str, Dict[str, Any]] = {}
    tree_samples: set = set()

    # Output file names derived from species and subtype. These are the same
    # names the files are given inside the Microreact bundle.
    prefix = make_prefix(args.species, args.subtype)
    logging.info("Using prefix '%s'", prefix)
    summary_out_file = f"{prefix}_summary.csv"
    matrix_out_file = f"{prefix}_dist.csv"
    tree_out_file = f"{prefix}.nwk"

    # ----------------------------
    # Alignment stats
    # ----------------------------
    stats_dict = list2map(load_csv(args.aln_stats, detect_sep(args.aln_stats)), key='name')
    if stats_dict:
        extend_dict(data, stats_dict)
    logging.info("Loaded alignment stats (%d rows)", len(stats_dict))

    # ----------------------------
    # SNP distance matrix (optional)
    # ----------------------------
    snp_dists = None
    if args.dist:
        snp_labels, snp_dists = load_dist(args.dist, matrix_out_file)
        extend_dict(data, find_links(snp_dists, snp_labels, args.strong_link, args.inter_link))
        # Preserve prior behavior: also write raw matrix as out.csv
        np.savetxt("out.csv", snp_dists, delimiter=",", fmt="%.6g")
        logging.info("Loaded distance matrix with %d samples", len(snp_labels))

    # ----------------------------
    # Tree (optional)
    # ----------------------------
    if args.tree:
        with open_text(args.tree) as f:
            tree = Phylo.read(f, 'newick')
        tree.root_at_midpoint()

        # Attempt to scale by reference length if present
        scale = 1.0
        try:
            vb = reference_row(stats_dict).get('length')
            if vb:
                scale = float(vb)
            else:
                logging.warning("Branch lengths not scaled!")
        except Exception:
            logging.warning("Branch lengths not scaled!")

        for clade in tree.find_clades():
            bl = getattr(clade, "branch_length", None)
            if bl is not None:
                try:
                    clade.branch_length = float(bl) * scale
                except Exception:
                    pass

        Phylo.write(tree, tree_out_file, "newick")

        labels, dist = patristic_distance_matrix(tree)
        tree_samples = set(labels)
        parts = partition_tree(dist, eps=args.partition_distance)
        extend_dict(data, {l: {'partition': p} for l, p in zip(labels, parts)})

        logging.info("Processed tree with %d tips and partitions", len(labels))

    # ----------------------------
    # Summary, filtered by alignment stats
    # ----------------------------
    original_fieldnames: List[str] = []
    samplesheet_ids: set = set()
    summary_rows = load_csv(args.summary, detect_sep(args.summary))
    if summary_rows:
        original_fieldnames = list(summary_rows[0].keys())
        samplesheet_ids = {r.get('sample', '') for r in summary_rows}

        filtered_rows = [r for r in summary_rows if r.get('sample', '') in stats_dict]
        logging.info("Filtered summary by alignment stats: %d -> %d", len(summary_rows), len(filtered_rows))

        summary_dict = list2map(filtered_rows, key='sample')
        if summary_dict:
            extend_dict(data, summary_dict)

    logging.info("Loaded and integrated summary file")

    # ----------------------------
    # Build final summary and write CSV (preserving requested column order)
    # ----------------------------
    summary_records: List[Dict[str, Any]] = []
    created_cols_order: List[str] = []   # tracks columns not in original_fieldnames
    created_seen: set = set()

    if data:
        include_only = tree_samples if tree_samples else None
        for sample_id, sample_data in data.items():
            if include_only and (sample_id not in include_only):
                continue
            rec = {'sample': sample_id}
            if is_reference_sample(sample_id):
                rec['status'] = ''
            elif sample_id in samplesheet_ids:
                rec['status'] = 'new'
            else:
                rec['status'] = 'old'
            for k, v in sample_data.items():
                if k in ('sample', 'status'):
                    continue
                if v != sample_id:  # avoid duplicating id if present as a value
                    rec[k] = ':'.join(map(str, v)) if isinstance(v, list) else v

                    # Track columns that are NEW (not from the original summary header)
                    if (k not in original_fieldnames) and (k not in created_seen):
                        created_cols_order.append(k)
                        created_seen.add(k)

            summary_records.append(rec)

    if summary_records:
        # Collect all keys observed
        all_keys = set()
        for rec in summary_records:
            all_keys.update(rec.keys())

        # 1) sample, then status
        fieldnames: List[str] = ['sample', 'status']

        # 2) newly created columns (in first-seen order)
        #    Keep only those that actually appear in the data
        fieldnames.extend([c for c in created_cols_order if c in all_keys and c not in ('sample', 'status')])

        # 3) remaining columns from the original summary, in their original order
        fieldnames.extend([fld for fld in original_fieldnames if fld not in ('sample', 'status') and fld in all_keys])

        with open(summary_out_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
            writer.writeheader()
            writer.writerows(summary_records)

        logging.info("Wrote %s with %d records and %d columns",
                    summary_out_file, len(summary_records), len(fieldnames))
    else:
        logging.warning("No records to write")

    # ----------------------------
    # Microreact packaging (optional)
    # ----------------------------
    if args.microreact:
        mr_json = load_json(args.microreact)
        mr_json.setdefault('meta', {})
        mr_json.setdefault('files', {})
        mr_json.setdefault('tables', {})
        mr_json['tables'].setdefault('table-1', {})  # ensure table exists

        # Attach files if present. The on-disk file names already use the
        # prefix, so they are reused directly as the Microreact attachment names.
        if _ensure_path(Path(summary_out_file), "Summary"):
            _attach_text_file(mr_json, "summary_file", Path(summary_out_file), summary_out_file)
        if _ensure_path(Path(tree_out_file), "Tree"):
            _attach_text_file(mr_json, "tree_file", Path(tree_out_file), tree_out_file)
        if _ensure_path(Path(matrix_out_file), "Distance matrix"):
            _attach_text_file(mr_json, "dist_file", Path(matrix_out_file), matrix_out_file)

        # Meta
        mr_json['meta']['name'] = make_meta_name(
            args.species, args.subtype, run_time.date().isoformat()
        )
        mr_json['meta']['timestamp'] = (
            run_time
            .isoformat(timespec="milliseconds")
            .replace("+00:00", "Z")
        )

        # Columns from written summary (if exists)
        if Path(summary_out_file).exists():
            with open(summary_out_file, "r", encoding="utf-8", newline="") as f:
                header = next(csv.reader(f), [])
            if header:
                mr_json['tables']['table-1']['columns'] = [{"field": h, "fixed": False} for h in header]

        out_path = f"{prefix}.microreact"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(mr_json, f, indent=2, ensure_ascii=False)
        logging.info("Wrote %s", out_path)

if __name__ == '__main__':
    main()