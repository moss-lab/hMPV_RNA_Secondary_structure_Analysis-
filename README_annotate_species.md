# annotate_species.py

**Stockholm species annotation**

Appends species names to the accessions in a Stockholm alignment, preserving all markup including the
#=GC SS_cons consensus structure line.

## Usage

```
python3 annotate_species.py alignment.sto --out alignment_named.sto
```

## Options

- `--root` — default `/lustre/hdd/LAS/wmoss-lab/ddabiri/CaCo_folder/MNV`
- `--motifs` — CSV/TSV of the motifs to process (needs Gene, Motif no., Coordinates). Defaults to the built-in list of the 29 covarying motifs; pass a file to override, or 'all' for every motif found; default `builtin`
- `--series` — restrict to this ScanFold run, e.g. NC_hMNV_cRNA
- `--out` — default `species_composition.csv`
- `--annotated-dir` — default `annotated_alignments`
- `--audit` — per-accession assignment and its source
- `--unresolved` — default `unresolved_accessions.txt`
- `--map` — TSV of accession <tab> species or organism text
- `--fetch` — look unknown accessions up at NCBI
- `--email` — contact address, required by NCBI
- `--api-key` — NCBI API key (raises the rate limit)
- `--cache` — default `ncbi_accessions.json`
- `--rename` — also append the species to each sequence name
- `--no-annotate` — build the table only, write no alignments

## Output

```
annotated Stockholm alignment with organism names prefixed to each accession
```

## Notes

Organism names come from the NCBI Entrez summary service; subgroup information is parsed from the sequence title.
