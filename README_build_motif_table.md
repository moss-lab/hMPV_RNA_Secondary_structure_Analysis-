# build_motif_table.py

**Motif summary table**

Builds the motif table from a ScanFold output folder. Covarying base-pair counts are read from the R-scape
.cacofold.power files rather than transcribed from summary output, after an earlier version of the table was
found to carry base-pair totals in the covarying column.

## Usage

```
python3 build_motif_table.py /path/to/scanfold_folder \
    --covariation cacofold-power --out motif_table.tsv
```

## Options

- `--root` — default `/lustre/hdd/LAS/wmoss-lab/ddabiri/CaCo_folder/MNV`
- `--out` — default `motif_table.tsv`
- `--markdown`
- `--audit`
- `--crosswalk`
- `--series` — keep only this ScanFold run, e.g. --series NC_hMNV_cRNA. Repeatable
- `--covariation` — where Covarying bp comes from (default: cacofold-power); one of `cacofold-power`, `given-power`, `cacofold-cov`, `given-cov`; default `cacofold-power`
- `--caco-pdf` — defaults to <root>/All.CaCo.pdf when present
- `--pdf-text` — use pre-extracted PDF text instead of the PDF
- `--dump-pdf-text`
- `--no-restrict` — keep every motif on disk, not only those in the PDF
- `--match-by` — one of `auto`, `coords`, `number`; default `auto`
- `--scanfold-table` — ScanFold summary table carrying MFE/z-score/ED
- `--compute-missing` — refold with ViennaRNA to fill missing MFE/ED
- `--covarying-out` — also write a table of only the motifs with covarying base pairs. With no value, the name is derived from --out (motif_table_covarying.tsv)
- `--min-covarying` — how many covarying pairs a motif needs to appear in that file (default 1); default `1`
- `--na` — what to write for a missing value (default NA; pass '' for an empty cell); default `NA`
- `--extra-columns` — also emit BP considered and Expected to covary
- `--explain` — show how this motif's covariation was derived, then exit
- `--debug`

## Output

```
motif_table.tsv   one row per motif with coordinates, MFE, z-score, ED and covariation
```

## Notes

Use --covariation cacofold-power. --series disambiguates folders holding more than one ScanFold run.
