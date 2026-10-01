#!/usr/bin/env python3

import urllib.request
import csv
import statistics
import os
import subprocess
import json
from pathlib import Path
import datetime

URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"

# Acceptable range = mean +/- Z_SCORE * stdev (2.58 ~ 99% of a normal distribution)
Z_SCORE = 2.58


def download(url, dest):
    """Download a file from URL to destination path."""
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        if getattr(r, "status", 200) != 200:
            raise RuntimeError(f"HTTP error {r.status} for {url}")
        while True:
            chunk = r.read(8192)
            if not chunk:
                break
            f.write(chunk)


def parse_float(val):
    """Parse a float value, returning None for invalid/missing values."""
    if not val or val.strip().lower() in {"na", "nan", "null", ""}:
        return None
    try:
        return float(val.strip())
    except ValueError:
        return None


def summarize(values, lower=None, upper=None):
    """
    Mean, standard deviation and acceptable range (mean +/- Z_SCORE * stdev)
    of a list of values. The range is clipped to [lower, upper] when given.
    Returns None if the stdev can't be calculated (fewer than 2 values).
    """
    if len(values) < 2:
        return None
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    low = mean - Z_SCORE * stdev
    high = mean + Z_SCORE * stdev
    if lower is not None:
        low = max(low, lower)
    if upper is not None:
        high = min(high, upper)
    return {"mean": mean, "stdev": stdev, "range": [low, high]}


def pull_taxids(taxids_file="taxids.txt", json_out="taxid_species.json"):
    """
    Query NCBI taxonomy database to get species-level information for taxids.
    Caches results in a JSON file to avoid re-querying.
    """
    taxids_path = Path(taxids_file)
    json_path = Path(json_out)

    # If JSON already exists, just load and return it
    if json_path.exists():
        with json_path.open() as f:
            return json.load(f)

    if not taxids_path.exists():
        raise FileNotFoundError(f"taxids file not found: {taxids_path}")

    # Build a bash script. Use the *Python* variable for the input path.
    # xtract output fields:
    #   0 taxid, 1 rank, 2 name, 3 lineage_species_taxid, 4 lineage_species_name
    cmd = rf"""
    set -euo pipefail

    epost -db taxonomy -input "{taxids_path}" \
      | efetch -format xml \
      | xtract \
          -pattern Taxon \
          -element TaxId Rank ScientificName \
          -block "LineageEx/Taxon" -if Rank -equals species -first TaxId ScientificName \
      | jq -Rn '
          [inputs
           | split("\t")
           | {{
               taxid: (.[0] | tonumber),
               species_taxid: (
                 if (.[1] // "") == "species"
                 then (.[0] | tonumber)
                 else (.[3] | tonumber?)
                 end
               ),
               species_name: (
                 if (.[1] // "") == "species"
                 then (.[2] // null)
                 else (.[4] // null)
                 end
               )
             }}
          ]'
    """

    # Run and write JSON
    with json_path.open("w") as outfile:
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            executable="/bin/bash",
            stdout=outfile,
            stderr=subprocess.PIPE,  # helpful error messages
            text=True,
        )

    # Load and return freshly created JSON
    with json_path.open() as f:
        return json.load(f)

def main():
    # Download or use existing eukaryotes data file
    src = 'eukaryotes.txt'
    if not os.path.exists(src):
        print(f"Downloading {URL} -> {src}")
        download(URL, src)
    else:
        print(f"Using existing {src}")

    # Parse the TSV file
    print("Parsing...")
    with open(src, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Parsed {len(rows)} records")

    # Build map of fungal taxids with their genome data
    taxid_map = {}
    for rec in rows:
        group = rec.get("Group")
        if group != "Fungi":
            continue
            
        taxid = rec.get("TaxID")
        if not taxid:
            continue
        
        taxid = str(taxid)
        name = rec.get('#Organism/Name')
        gc = parse_float(rec.get("GC%"))
        length_mb = parse_float(rec.get("Size (Mb)"))

        if taxid not in taxid_map:
            taxid_map[taxid] = {
                "name": [],
                "length": [],
                "gc": []
            }

        if name:
            taxid_map[taxid]['name'].append(name)
        if length_mb is not None:
            taxid_map[taxid]['length'].append(length_mb)
        if gc is not None:
            taxid_map[taxid]['gc'].append(gc)

    # Write taxids to file for NCBI query
    with open('taxids.txt', 'w') as f:
        f.write('\n'.join(taxid_map.keys()))

    # Get species-level taxonomy information
    print("Querying NCBI taxonomy...")
    species_map = pull_taxids()

    # Roll up data to species level
    final_map = {}
    for rec in species_map:
        taxid = str(rec.get('taxid'))
        species_taxid = rec.get('species_taxid')
        species_name = rec.get('species_name')

        # Skip records without complete information
        if not taxid or not species_taxid or not species_name:
            continue
        
        # Skip taxids we don't have genome data for
        if taxid not in taxid_map:
            continue

        # Initialize species entry if needed
        if species_taxid not in final_map:
            final_map[species_taxid] = {
                'species_name': species_name,
                'taxids': [str(species_taxid)],
                'names': [species_name],
                'length': [],
                'gc': []
            }
        
        # Add data from this taxid to the species
        final_map[species_taxid]['taxids'].append(taxid)
        final_map[species_taxid]['names'].extend(taxid_map[taxid]['name'])
        final_map[species_taxid]['length'].extend(taxid_map[taxid]['length'])
        final_map[species_taxid]['gc'].extend(taxid_map[taxid]['gc'])

    # Calculate statistics for each species
    print("Calculating statistics...")
    species_stats = []
    skipped = 0
    for species_taxid, data in final_map.items():
        # Filter out None values
        lengths = [float(l) * 10**6 for l in data['length'] if l is not None]
        gcs = [float(g) for g in data['gc'] if g is not None]
        
        # Mean, stdev and range; None if fewer than 2 values
        length = summarize(lengths, lower=0)
        gc = summarize(gcs, lower=0, upper=100)

        # Keep only species where both ranges could be determined
        if length is None or gc is None:
            skipped += 1
            continue

        # Get unique names and taxids
        unique_names = sorted(set([n for n in data['names'] if n]))
        unique_taxids = sorted(set(data['taxids']))

        # Calculate stats
        stats = {
            'species_name': data['species_name'],
            'species_taxid': species_taxid,
            'taxids': unique_taxids,
            'names': unique_names,
            'n': len(lengths),
            'z_score': Z_SCORE,
            'length_mean': length['mean'],
            'length_stdev': length['stdev'],
            'length_range': length['range'],
            'gc_mean': gc['mean'],
            'gc_stdev': gc['stdev'],
            'gc_range': gc['range']
        }
        species_stats.append(stats)
    
    # Sort by species name
    species_stats.sort(key=lambda x: x['species_name'])
    
    # Write output JSON
    timestamp = datetime.date.today().isoformat()
    output_file = f'{timestamp}_ncbi-fungal-sp.json'
    with open(output_file, 'w') as f:
        json.dump(species_stats, f, indent=2)
    
    print(f"Wrote statistics for {len(species_stats)} species to {output_file}")
    print(f"Skipped {skipped} species with fewer than 2 genomes with length and GC data")


if __name__ == "__main__":
    main()