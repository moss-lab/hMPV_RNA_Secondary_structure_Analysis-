# make_titv_table.py

**Transition/transversion table (earlier version)**

Earlier transition/transversion table builder, kept for reference. Superseded by titv_table.py, which adds
significance testing, confidence intervals and codon-position restriction.

## Usage

```
python3 make_titv_table.py <alignments> --out titv_old.tsv
```

## Options

- `summary` — titv_summary.csv from titv_analysis.py
- `--metadata` — CSV/TSV of motif gene/coordinates/MFE/z-score/ED
- `--out` — default `table4_titv.csv`
- `--markdown` — also write an aligned markdown table
- `--q-cut` — default `0.05`
- `--min-subs` — below this a row is flagged low-n (default 20); default `20`
- `--clean-only` — keep only rows with no flags
- `--sort` — one of `motif`, `q`, `or`; default `motif`
- `--na` — default `NA`

## Output

```
titv_old.tsv   stem and loop ratios per structure
```

## Notes

Use titv_table.py instead.
