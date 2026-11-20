# resolve_accessions.py

**Accession to organism lookup**

Resolves NCBI nucleotide accessions to organism names. Standard library only, no Biopython required.

## Usage

```
python3 resolve_accessions.py accessions.txt --out acc_species.tsv
```

## Options

- `accessions` — file with one accession per line
- `--email` — contact address, required by NCBI; required
- `--api-key` — NCBI API key, if you have one
- `--out` — default `accession_species.tsv`
- `--cache` — default `ncbi_records.json`

## Output

```
acc_species.tsv   accession and organism name
```

## Notes

Queries NCBI in batches and caches results, so it can be rerun without repeating lookups.
