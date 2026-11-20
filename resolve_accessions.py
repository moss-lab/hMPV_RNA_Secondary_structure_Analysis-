#!/usr/bin/env python3
"""Resolve NCBI nucleotide accessions to organism names."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

ACCESSION_RE = re.compile(r"\b([A-Z]{1,2}_?\d{5,8}(?:\.\d+)?)\b")

# Bovine before the generic RSV patterns; subgroup A/B before unsubtyped.
PATTERNS = [
    ("Bovine_RSV", re.compile(r"bovine|\bBRSV\b", re.IGNORECASE)),
    ("hMPV", re.compile(r"metapneumovirus|\bhMPV\b", re.IGNORECASE)),
    ("RSV-A", re.compile(
        r"(?:syncytial\s+virus|orthopneumovirus|\bRSV\b|\bHRSV\b)"
        r"[^A-Za-z0-9]{0,25}(?:subgroup\s+|subtype\s+|type\s+|group\s+)?A\b"
        r"|/A/|\bRSVA\b|\bsubgroup\s+A\b|\btype\s+A\b", re.IGNORECASE)),
    ("RSV-B", re.compile(
        r"(?:syncytial\s+virus|orthopneumovirus|\bRSV\b|\bHRSV\b)"
        r"[^A-Za-z0-9]{0,25}(?:subgroup\s+|subtype\s+|type\s+|group\s+)?B\b"
        r"|/B/|\bRSVB\b|\bsubgroup\s+B\b|\btype\s+B\b", re.IGNORECASE)),
    ("RSV-unassigned", re.compile(
        r"respiratory\s+syncytial\s+virus|orthopneumovirus|\bRSV\b",
        re.IGNORECASE)),
]

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def classify(text: str) -> str:
    for label, pattern in PATTERNS:
        if pattern.search(text or ""):
            return label
    return "UNCLASSIFIED"


def read_accessions(path: str) -> list[str]:
    seen, out = set(), []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = ACCESSION_RE.search(line)
            acc = (m.group(1) if m else line).split(".")[0]
            if acc not in seen:
                seen.add(acc)
                out.append(acc)
    return out


def fetch(accessions, email, api_key=None, batch=180, cache_path=None):
    cache = {}
    if cache_path and os.path.isfile(cache_path):
        try:
            cache = json.load(open(cache_path))
            print(f"Reusing {len(cache)} cached records")
        except (OSError, ValueError):
            cache = {}

    todo = [a for a in accessions if a not in cache]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        params = {"db": "nuccore", "id": ",".join(chunk), "retmode": "json",
                  "tool": "resolve-accessions", "email": email}
        if api_key:
            params["api_key"] = api_key
        print(f"  requesting {i + 1}-{i + len(chunk)} of {len(todo)}...",
              flush=True)
        for attempt in (1, 2, 3):
            try:
                with urllib.request.urlopen(
                        f"{ESUMMARY}?{urllib.parse.urlencode(params)}",
                        timeout=90) as resp:
                    data = json.load(resp)
                break
            except Exception as exc:
                if attempt == 3:
                    print(f"  ! giving up on this batch: {exc}",
                          file=sys.stderr)
                    data = {}
                else:
                    print(f"  ! {exc}; retrying in {attempt * 5}s",
                          file=sys.stderr)
                    time.sleep(attempt * 5)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            acc = (rec.get("accessionversion") or rec.get("caption")
                   or "").split(".")[0]
            text = " ".join(filter(None, [rec.get("organism", ""),
                                          rec.get("title", ""),
                                          rec.get("strain", "")])).strip()
            if acc:
                cache[acc] = text
        if cache_path:
            try:
                json.dump(cache, open(cache_path, "w"), indent=1,
                          sort_keys=True)
            except OSError:
                pass
        time.sleep(0.4 if api_key else 1.0)   # NCBI rate limit
    return cache


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accessions", help="file with one accession per line")
    ap.add_argument("--email", required=True,
                    help="contact address, required by NCBI")
    ap.add_argument("--api-key", help="NCBI API key, if you have one")
    ap.add_argument("--out", default="accession_species.tsv")
    ap.add_argument("--cache", default="ncbi_records.json")
    args = ap.parse_args()

    accs = read_accessions(args.accessions)
    print(f"{len(accs)} distinct accessions to resolve")

    records = fetch(accs, args.email, args.api_key, cache_path=args.cache)

    counts = {}
    with open(args.out, "w") as fh:
        fh.write("# accession\tspecies\torganism text from NCBI\n")
        for acc in accs:
            text = records.get(acc, "")
            label = classify(text) if text else "UNCLASSIFIED"
            counts[label] = counts.get(label, 0) + 1
            fh.write(f"{acc}\t{label}\t{text or 'NOT FOUND'}\n")

    print(f"\nWrote {args.out}")
    for label in sorted(counts, key=lambda k: -counts[k]):
        print(f"    {label:16s} {counts[label]}")
    if counts.get("UNCLASSIFIED"):
        print(f"\n{counts['UNCLASSIFIED']} accessions did not match any "
              f"pattern. Look at their organism text in the third column: "
              f"either the record is something unexpected, or the pattern "
              f"needs widening. They stay UNCLASSIFIED until you decide.")


if __name__ == "__main__":
    main()
