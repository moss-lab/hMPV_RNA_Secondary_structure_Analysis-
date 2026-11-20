# export_motifs.py

**Motif structure export**

Writes motif_N.dbn and motif_N.sto for a list of motifs. Converts WUSS notation to plain dot-bracket and
rebalances brackets left orphaned when gapped columns are removed.

## Usage

```
python3 export_motifs.py --folder /path/to/scanfold_folder --out exported_motifs
```

## Options

- `--root` — folder holding the alignments and structures (default: the current directory); default `.`
- `--out` — output directory (default exported_motifs); default `exported_motifs`
- `--series` — restrict to one ScanFold run, e.g. NC_hMNV_cRNA
- `--motifs` — CSV/TSV or plain list of motif numbers to export, overriding the built-in default set
- `--all-motifs` — export every motif found, ignoring the built-in default set
- `--prefix` — output filename prefix (default 'motif_', giving motif_13.dbn and motif_13.sto); default `motif_`
- `--sto-format` — fasta (default) writes the >accession layout into .sto; stockholm writes real Stockholm; both writes .sto and .stockholm; one of `fasta`, `stockholm`, `both`; default `fasta`
- `--use-cacofold` — prefer the .cacofold alignment and structure
- `--ss-from-cov` — if no .dbn and no SS_cons, rebuild the structure from the R-scape pair list
- `--keep-pseudoknots` — keep crossing pairs as [] brackets; by default they are written as unpaired, since plain dot-bracket cannot represent them
- `--ref` — substring of the sequence the .dbn describes, used when mapping it onto alignment columns
- `--degap-to` — write the structure degapped to this sequence (substring match) instead of to alignment columns
- `--dbn-with-sequence` — include the reference sequence line in the .dbn
- `--dbn-header` — write a '>motif_N' header line into the .dbn
- `--names` — sequence header style (default accession); one of `accession`, `species`, `both`, `raw`; default `accession`
- `--name-map` — two-column accession->label file for --names
- `--no-version` — strip the .1 version suffix from accessions
- `--wrap` — wrap sequences at N characters (0 = one line)
- `--min-seqs` — skip motifs with fewer sequences than this; default `1`
- `--manifest` — per-motif record of what was written; default `export_manifest.tsv`
- `--dry-run` — report what would be written, create nothing

## Output

```
exported_motifs/motif_N.dbn   dot-bracket structure
exported_motifs/motif_N.sto   alignment, aligned FASTA or Stockholm
```

## Notes

--sto-format selects fasta, stockholm or both. Sequence names containing spaces are sanitised, since NCFBI
FASTA descriptions otherwise parse as sequence and corrupt the alignment silently.
