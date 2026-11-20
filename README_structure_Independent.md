# structure_agnostic.py

**Whole-genome structure-agnostic test (superseded)**

Whole-genome version of the structure-agnostic covariation test. Superseded by pergene_agnostic.py once it
became clear that the differing gene order between Metapneumovirus and Orthopneumovirus makes whole-genome
alignment across the dataset meaningless.

## Usage

```
python3 structure_agnostic.py --fasta <genomes> --motifs <table> --out <dir>
```

## Options

- `--genomes` — FASTA of related genomes, including the reference; required
- `--reference` — accession or name substring of the genome whose coordinates the motif table uses, e.g. NC_039199; required
- `--motifs` — table with motif numbers and coordinates; required
- `--structures` — folder of per-motif .dbn files. Without it the regions are extracted but no structure is tested
- `--out` — default `agnostic`
- `--alignment` — use this existing genome alignment instead of running MAFFT
- `--mafft-preset` — default `auto`
- `--mafft-extra` — default `--quiet --thread -1`
- `--max-genomes` — subsample to this many genomes (the reference is always kept). Use for a trial run
- `--flank` — extra columns either side of each motif
- `--min-residues` — drop a sequence from a motif region if it has fewer than this many non-gap residues there; default `10`
- `--rscape-bin` — path to the R-scape executable if it is not simply 'R-scape' on PATH; default `R-scape`
- `--mafft-bin` — path to the MAFFT executable; default `mafft`
- `--rscape-args` — R-scape flags. Use the SAME settings as your original analysis, otherwise a difference in the covariation counts could come from the settings rather than from the alignment. Drop --nofigures to have R-scape generate R2R figures; default `-s --nofigures`
- `--compare` — earlier results table to place beside the new counts (needs Motif and a covariation column)
- `--dry-run`

## Output

```
same table set as pergene_agnostic.py, computed on a single whole-genome alignment
```

## Notes

Not used for the published analysis. Use pergene_agnostic.py.
