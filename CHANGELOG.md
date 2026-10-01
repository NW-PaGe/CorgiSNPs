# NW-PaGe/corgisnps: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of NW-PaGe/corgisnps, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Variant calling and phylogenetic parameters can be set per species or per subtype in the reference database manifest. Precedence is subtype > species > run-level parameter; overrides are validated on load and listed in the log.
- Automated QC checks de novo assembly length and GC content against a min / max range per species instead of a z-score. Ranges can be set per species in the reference manifest (`length_range`, `gc_range`) and are otherwise taken from the NCBI statistics (`make_ncbi_stats.py` now pre-computes them). Summary columns `denovo_length_range`, `denovo_gc_range` and `qc_range_source` replace `denovo_length_z` and `denovo_gc_z`.
- Species aliases in the reference manifest are used to find a species in the NCBI statistics (e.g. a samplesheet species of `Candida auris`).

### `Fixed`

### `Dependencies`

### `Deprecated`

- `--max_z_score_qc` has been removed; assembly QC now uses the species' length and GC ranges.
