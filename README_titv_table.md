# titv_table.py

**Transition/transversion analysis**

Compares transition and transversion counts between paired and unpaired positions of each structure.
Sequences are Henikoff-weighted, counts are pooled across columns before any ratio is formed, significance
is assessed by two-sided Fisher exact test, and P values are corrected across structures by
Benjamini-Hochberg. The neutral expectation for either ratio is 0.5, not 1.

## Usage

```
python3 titv_table.py pergene/motif_*/*_pergene.sto --out titv --weight --markdown
```

## Options

- `inputs` — alignment files (.sto/.fa/...); .dbn files may be mixed in and will be paired by name
- `--dbn-dir` — directory of .dbn files
- `--dbn` — one .dbn applied to a single alignment
- `--out` — output directory (default titv_results); default `titv_results`
- `--table` — final table filename inside --out; default `titv_table.csv`
- `--markdown` — also write the table as aligned markdown
- `--per-structure` — also write a per-mutation-type CSV per structure
- `--restrict` — all columns (default), third codon positions, or fourfold-degenerate sites. Use one of the latter for coding regions to remove protein-level constraint; one of `all`, `third`, `fourfold`; default `all`
- `--frame` — one of `0`, `1`, `2`
- `--frames` — CSV with columns name,frame
- `--ref` — substring of the sequence the structure describes, for mapping onto columns
- `--weight` — Henikoff sequence weighting
- `--boot` — default `5000`
- `--seed`
- `--strip` — extra regex stripped from filenames when pairing an alignment with its structure. Repeatable; default `[]`
- `--metadata` — CSV/TSV joined into the final table
- `--meta-key` — column in --metadata holding the join key; default `Motif`
- `--meta-match` — auto matches on the trailing number in each name (motif_12 <-> 12); name matches the whole string; one of `auto`, `name`; default `auto`
- `--meta-columns` — comma-separated metadata columns to place at the front of the table (default: all of them)
- `--columns` — comma-separated list overriding the output columns entirely. Any field computed by the script is available; --list-fields prints them
- `--list-fields` — print every available column name and exit
- `--q-cut` — default `0.05`
- `--min-subs` — rows with fewer substitutions are flagged low-n; default `20`
- `--sort` — one of `input`, `key`, `q`, `or`; default `key`

## Output

```
titv.tsv   structure, n_seq, n_subs, stem and loop ti/tv, odds ratio, CI, P, q, direction
```

## Notes

--restrict third and --restrict fourfold limit the comparison to third or fourfold-degenerate codon positions.
The fourfold analysis has not been run and is needed: in coding sequence the stem/loop asymmetry may reflect
codon position rather than base pairing.
Odds ratios computed from a class containing no transversions are unstable and are flagged rather than reported.
