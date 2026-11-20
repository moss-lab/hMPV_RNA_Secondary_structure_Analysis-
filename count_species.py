#!/usr/bin/env python3
"""Count sequences per species in one alignment or a folder of alignments."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

ALN_EXT = (".sto", ".stockholm", ".stk", ".fa", ".fasta", ".fna", ".afa",
           ".aln", ".afasta")
ACCESSION_RE = re.compile(r"\b([A-Z]{1,2}_?\d{5,8}(?:\.\d+)?)\b")
MOTIF_RE = re.compile(r"motif[_-](\d+)(?=[._\-]|$)", re.IGNORECASE)

# Order matters. Bovine is tested before the generic RSV patterns, and the
# subgroup patterns before unsubtyped RSV, so a name reading
# "bovine respiratory syncytial virus" is never counted as HRSV.
RULES = [
    ("BRSV", re.compile(
        r"bovine|\bBRSV\b|\bbRSV\b|Bovine_RSV|bovine[_\s-]*orthopneumo",
        re.IGNORECASE)),
    ("hMPV", re.compile(
        r"metapneumovirus|\bhMPV\b|\bHMPV\b|\bAMPV\b|\bMPV\b",
        re.IGNORECASE)),
    ("HRSV-A", re.compile(
        r"(?:hRSV|HRSV|RSV|syncytial[_\s]*virus|orthopneumovirus)"
        r"[^A-Za-z0-9]{0,25}(?:subgroup|subtype|type|group)?[^A-Za-z0-9]{0,6}A"
        r"(?![A-Za-z])"
        r"|\bRSVA\b|\bRSV[_-]A\b|/A/|\bsubgroup[_\s]*A\b",
        re.IGNORECASE)),
    ("HRSV-B", re.compile(
        r"(?:hRSV|HRSV|RSV|syncytial[_\s]*virus|orthopneumovirus)"
        r"[^A-Za-z0-9]{0,25}(?:subgroup|subtype|type|group)?[^A-Za-z0-9]{0,6}B"
        r"(?![A-Za-z])"
        r"|\bRSVB\b|\bRSV[_-]B\b|/B/|\bsubgroup[_\s]*B\b",
        re.IGNORECASE)),
    ("RSV-?", re.compile(
        r"respiratory[_\s]*syncytial|orthopneumovirus|\bRSV\b|\bHRSV\b",
        re.IGNORECASE)),
]

COLUMNS = ["hMPV", "HRSV-A", "HRSV-B", "BRSV", "RSV-?", "Unknown"]


def classify(text):
    for label, pattern in RULES:
        if pattern.search(text or ""):
            return label
    return None


def sniff(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            t = line.strip()
            if not t:
                continue
            if t.startswith("# STOCKHOLM"):
                return "stockholm"
            if t.startswith(">"):
                return "fasta"
            if t.startswith("#"):
                continue
            return "stockholm" if len(t.rsplit(None, 1)) == 2 else "fasta"
    return "fasta"


def sequence_names(path):
    """Distinct sequence names, in order. Interleaved Stockholm blocks repeat."""
    names, seen = [], set()
    fmt = sniff(path)
    gs = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            s = line.rstrip("\n").strip()
            if not s or s == "//":
                continue
            if fmt == "fasta":
                if s.startswith(">"):
                    nm = s[1:].strip()
                    if nm not in seen:
                        seen.add(nm)
                        names.append(nm)
                continue
            # Stockholm: capture #=GS DE descriptions, skip other markup
            if s.startswith("#=GS"):
                parts = s.split(None, 3)
                if len(parts) >= 4 and parts[2].upper() == "DE":
                    gs[parts[1]] = parts[3].strip()
                continue
            if s.startswith("#"):
                continue
            parts = s.rsplit(None, 1)
            if len(parts) != 2:
                continue
            nm = parts[0].strip()
            if nm not in seen:
                seen.add(nm)
                names.append(nm)
    return names, gs


def load_map(path):
    out = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]", line, maxsplit=1)
            if len(parts) == 2:
                out[parts[0].strip().split(".")[0]] = parts[1].strip()
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None,
                    help="a folder of alignments, or a single alignment "
                         "file. Same as --root.")
    ap.add_argument("--root", default=None,
                    help="folder of alignments (default: current directory)")
    ap.add_argument("--no-recursive", action="store_true")
    ap.add_argument("--out", default="species_counts.csv",
                    help="per-file counts (.csv or .tsv)")
    ap.add_argument("--audit", help="write every sequence name and its call")
    ap.add_argument("--map", dest="map_file",
                    help="accession<tab>species, for names lacking organism "
                         "text")
    ap.add_argument("--list-unknown", action="store_true",
                    help="print unresolved names so the pattern can be seen")
    ap.add_argument("--percent", action="store_true",
                    help="add a percentage column per category")
    ap.add_argument("--min-seqs", type=int, default=1)
    args = ap.parse_args()

    if args.folder and args.root and args.folder != args.root:
        sys.exit(f"Two folders given: {args.folder!r} and {args.root!r}. "
                 f"Pass one.")
    args.root = args.root or args.folder or "."
    single_file = os.path.isfile(args.root)
    if not single_file and not os.path.isdir(args.root):
        sys.exit(f"No such file or directory: {args.root}")

    files = []
    if single_file:
        # A single alignment is analysed whatever its extension: the
        # format is decided by content, not by the filename.
        files = [args.root]
    else:
        walker = (os.walk(args.root) if not args.no_recursive else
                  [(args.root, [], os.listdir(args.root))])
        for d, _sub, names in walker:
            for n in names:
                if n.lower().endswith(ALN_EXT):
                    files.append(os.path.join(d, n))
    files.sort(key=lambda p: (
        int(MOTIF_RE.search(os.path.basename(p)).group(1))
        if MOTIF_RE.search(os.path.basename(p)) else 10**9,
        os.path.basename(p)))
    if not files:
        sys.exit(f"No alignment files found under {args.root}. Recognised "
                 f"extensions: {', '.join(ALN_EXT)}. Pass a single file "
                 f"directly to analyse it whatever its name.")
    print(f"{len(files)} alignment file(s)" +
          ("" if single_file else f" under {args.root}"))

    name_map = load_map(args.map_file) if args.map_file else {}
    if args.map_file:
        print(f"Loaded {len(name_map)} accession mappings")

    rows, audit, totals = [], [], Counter()
    unknown_names, per_label_files = Counter(), defaultdict(set)

    for path in files:
        base = os.path.basename(path)
        try:
            names, gs = sequence_names(path)
        except OSError as exc:
            print(f"  ! cannot read {base}: {exc}", file=sys.stderr)
            continue
        if len(names) < args.min_seqs:
            continue

        counts = Counter()
        for nm in names:
            # the name itself, then its #=GS DE line, then the map
            label = classify(nm) or classify(gs.get(nm, ""))
            source = "name" if classify(nm) else (
                "#=GS DE" if classify(gs.get(nm, "")) else None)
            if label is None:
                m = ACCESSION_RE.search(nm)
                if m and name_map:
                    acc = m.group(1).split(".")[0]
                    label = (name_map[acc] if name_map.get(acc) in COLUMNS
                             else classify(name_map.get(acc, "")))
                    if label:
                        source = "--map"
            if label is None:
                label, source = "Unknown", "unresolved"
                unknown_names[nm] += 1
            counts[label] += 1
            per_label_files[label].add(base)
            audit.append({"file": base, "sequence": nm,
                          "species": label, "source": source})

        mm = MOTIF_RE.search(base)
        row = {"File": base, "Motif": int(mm.group(1)) if mm else ""}
        for c in COLUMNS:
            row[c] = counts.get(c, 0)
        row["Total"] = len(names)
        if args.percent:
            for c in COLUMNS:
                row[c + " %"] = (f"{100*counts.get(c,0)/len(names):.1f}"
                                 if names else "0.0")
        rows.append(row)
        totals.update(counts)

    if not rows:
        sys.exit("Nothing to count.")

    cols = (["File", "Motif"] + COLUMNS + ["Total"]
            + ([c + " %" for c in COLUMNS] if args.percent else []))
    delim = "\t" if args.out.lower().endswith((".tsv", ".txt")) else ","
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter=delim,
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote per-file counts to {args.out}")

    if args.audit:
        with open(args.audit, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["file", "sequence", "species",
                                               "source"],
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(audit)
        print(f"Wrote per-sequence calls to {args.audit}")

    grand = sum(totals.values())
    print(f"\n{'Category':12}{'Sequences':>11}{'Percent':>9}"
          f"{'Files with it':>15}")
    print("-" * 47)
    for c in COLUMNS:
        n = totals.get(c, 0)
        if n or c in ("hMPV", "HRSV-A", "HRSV-B", "BRSV"):
            print(f"{c:12}{n:>11}{100*n/grand if grand else 0:>8.1f}%"
                  f"{len(per_label_files[c]):>15}")
    print("-" * 47)
    print(f"{'TOTAL':12}{grand:>11}{100.0:>8.1f}%{len(rows):>15}")

    rsv = sum(totals.get(c, 0) for c in ("HRSV-A", "HRSV-B", "BRSV", "RSV-?"))
    print(f"\nhMPV {totals.get('hMPV',0)} vs RSV-group {rsv} "
          f"({100*rsv/grand if grand else 0:.1f}% of all sequences)")

    depths = [r["Total"] for r in rows]
    print(f"Alignment depth: min {min(depths)}, median "
          f"{sorted(depths)[len(depths)//2]}, max {max(depths)}")
    no_hmpv = [r for r in rows if r["hMPV"] == 0]
    if no_hmpv:
        print(f"\n  ! {len(no_hmpv)} alignment(s) contain NO hMPV sequence: "
              f"{', '.join(str(r['Motif'] or r['File']) for r in no_hmpv[:10])}",
              file=sys.stderr)
        print(f"    An hMPV structure with no hMPV sequence in its own "
              f"alignment needs checking.", file=sys.stderr)
    only_hmpv = [r for r in rows if r["Total"] == r["hMPV"]]
    if only_hmpv:
        print(f"  ! {len(only_hmpv)} alignment(s) are hMPV-only, so they "
              f"carry no cross-species covariation signal: "
              f"{', '.join(str(r['Motif'] or r['File']) for r in only_hmpv[:10])}",
              file=sys.stderr)
    if totals.get("RSV-?"):
        print(f"  ! {totals['RSV-?']} sequence(s) are human RSV with no "
              f"subgroup in the name. They are counted under 'RSV-?' and "
              f"are not assigned to A or B.", file=sys.stderr)
    if totals.get("Unknown"):
        print(f"  ! {totals['Unknown']} sequence(s) ({len(unknown_names)} "
              f"distinct names) could not be classified.", file=sys.stderr)
        if args.list_unknown:
            for nm, n in unknown_names.most_common():
                print(f"      {nm}   x{n}", file=sys.stderr)
        else:
            for nm, n in unknown_names.most_common(5):
                print(f"      e.g. {nm}   x{n}", file=sys.stderr)
            print(f"    Use --list-unknown for all of them, or --map to "
                  f"supply accession<tab>species.", file=sys.stderr)


if __name__ == "__main__":
    main()
