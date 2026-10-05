# myco_corgi_snpdists

Tools for converting MycoSNP VCFs to CorgiSNPs-style wide CSVs and comparing the two.

Requires Python 3.8+ and `vcfTools.py` (same directory) for the VCF converter. Use `environment.yml` to create the conda environment.

## Using Wrapper Script

### `run_mycosnp_vs_corgisnps.py`

Runs the VCF→wide CSV conversion, then the comparison. Accepts gzipped or plain reference FASTA and MycoSNP finalfiltered.vcf. Sets `--max_amb_samples` to 10% of the VCF sample count (floored, minimum 1).

```bash
python run_mycosnp_vs_corgisnps.py \
  --reference ref.fasta.gz \
  --vcf mycosnp-finalfiltered.vcf.gz \
  --corgisnps full.csv
```

Defaults write `mycosnp_full.csv`, `snp_comparison_summary.csv`, and `snp_comparison_details.csv`. Override paths as needed:

```bash
python run_mycosnp_vs_corgisnps.py \
  --reference ref.fasta \
  --vcf mycosnp-finalfiltered.vcf \
  --corgisnps full.csv \
  --mycosnp-csv mycosnp_full.csv \
  --output snp_comparison_summary.csv \
  --details-output snp_comparison_details.csv \
  --min_depth 10
```

## Individual scripts

### `vcfSnpsToFasta_genomewide_fullcsv.py`

Convert MycoSNP finalfiltered.vcf to a genome-wide wide CSV (`CHROM,POS,FILTER,REF,<samples>...`). Writes one row per reference position to **stdout**. Low-depth, filtered, or ambiguous calls become `N`. Plain (uncompressed) inputs only. Pass the name of the fourth column of CorgiSNPs full.csv as `--reference-column-name`, so they match for the next step.

```bash
python vcfSnpsToFasta_genomewide_fullcsv.py mycosnp-finalfiltered.vcf \
  --reference ref.fasta \
  --max_amb_samples 3 \
  --reference-column-name REF_NAME_FROM_CORGISNPS_FULL_CSV \
  > mycosnp_full.csv
```

### `compare_mycosnp_corgisnps.py`

Compare sample-wise SNP calls between a MycoSNP wide CSV and a CorgiSNPs `full.csv`. Both files must share coordinates, reference-column name, and sample names. Writes a per-sample summary CSV; optionally streams per-difference details.

```bash
python compare_mycosnp_corgisnps.py \
  --mycosnp mycosnp_full.csv \
  --corgisnps full.csv \
  --output snp_comparison_summary.csv \
  --details-output snp_comparison_details.csv \
  --progress-every 1000000
```
