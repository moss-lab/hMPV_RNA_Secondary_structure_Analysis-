# count_species.py

**Species composition**

Counts sequences per species in a single alignment or across a folder of alignments.

## Usage

```
python3 count_species.py Group_1.fa --out species_count.csv
```

## Options

- `folder` — a folder of alignments, or a single alignment file. Same as --root
- `--root` — folder of alignments (default: current directory)
- `--no-recursive`
- `--out` — per-file counts (.csv or .tsv); default `species_counts.csv`
- `--audit` — write every sequence name and its call
- `--map` — accession<tab>species, for names lacking organism text
- `--list-unknown` — print unresolved names so the pattern can be seen
- `--percent` — add a percentage column per category
- `--min-seqs` — default `1`

## Output

```
species_count.csv   counts of hMPV, RSV-A, RSV-B, bovine RSV, unassigned and unclassified
```

## Notes

Accepts either a file or a directory as the positional argument.
