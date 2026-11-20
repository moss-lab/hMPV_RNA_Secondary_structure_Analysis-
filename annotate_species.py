#!/usr/bin/env python3
"""Append species names to the accessions in a Stockholm alignment."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from typing import Optional

# Motif identification (same anchoring rules as build_motif_table.py)

MOTIF_RE = re.compile(r"motif[_-](\d+)(?=[._\-]|$)", re.IGNORECASE)
COORD_RE = re.compile(r"coordinates[:_](\d+)[-_](\d+)", re.IGNORECASE)


def parse_name(basename: str) -> tuple[Optional[str], Optional[int]]:
    m = MOTIF_RE.search(basename)
    if not m:
        return None, None
    series = basename[:m.start()].rstrip("_-.") or "(unnamed)"
    return series, int(m.group(1))


# Accessions and species

# GenBank/RefSeq style: NC_039199.1, KX027331.1, AF209337, JX198021.1
ACCESSION_RE = re.compile(
    r"\b([A-Z]{1,2}_?\d{5,8}(?:\.\d+)?)\b")

# Columns in the output, in order. "Unclassified" is appended only when
# something lands in it.
SPECIES_COLUMNS = ["hMPV", "RSV-A", "RSV-B", "Bovine_RSV"]

# Patterns applied to any organism text we have (sequence name, map file,
# or NCBI title). Order matters: bovine is tested before the generic RSV
# patterns, and subgroup A/B before unsubtyped RSV.
ORGANISM_PATTERNS = [
    ("Bovine_RSV", re.compile(
        r"bovine|\bBRSV\b|\bbRSV\b|respiratory\s+syncytial\s+virus\s+"
        r"type\s+[AB]?\s*bovine", re.IGNORECASE)),
    ("hMPV", re.compile(
        r"metapneumovirus|\bhMPV\b|\bHMPV\b|\bAMPV\b|avian\s+metapneumo",
        re.IGNORECASE)),
    ("RSV-A", re.compile(
        r"(?:syncytial\s+virus|orthopneumovirus|\bRSV\b|\bHRSV\b)"
        r"[^A-Za-z0-9]{0,20}(?:subgroup\s+|subtype\s+|type\s+|group\s+)?A\b"
        r"|/A/|\bRSVA\b|\bsubgroup\s+A\b", re.IGNORECASE)),
    ("RSV-B", re.compile(
        r"(?:syncytial\s+virus|orthopneumovirus|\bRSV\b|\bHRSV\b)"
        r"[^A-Za-z0-9]{0,20}(?:subgroup\s+|subtype\s+|type\s+|group\s+)?B\b"
        r"|/B/|\bRSVB\b|\bsubgroup\s+B\b", re.IGNORECASE)),
    # Unsubtyped human RSV: recorded distinctly rather than guessed into A.
    ("RSV-unassigned", re.compile(
        r"respiratory\s+syncytial\s+virus|orthopneumovirus|\bRSV\b|\bHRSV\b",
        re.IGNORECASE)),
]

# Reference genomes, so the common cases resolve with no network at all.
BUILTIN_ACCESSIONS = {
    "NC_039199": "hMPV",       # Human metapneumovirus
    "NC_004148": "hMPV",       # Human metapneumovirus (older RefSeq)
    "NC_038235": "RSV-A",      # Human orthopneumovirus, subgroup A
    "NC_001781": "RSV-B",      # Human orthopneumovirus, subgroup B
    "NC_001989": "Bovine_RSV", # Bovine orthopneumovirus
    "NC_038272": "Bovine_RSV",
}


def classify_organism(text: str) -> Optional[str]:
    """Species bucket from free organism text, or None."""
    if not text:
        return None
    for label, pattern in ORGANISM_PATTERNS:
        if pattern.search(text):
            return label
    return None


def bare_accession(acc: str) -> str:
    return acc.split(".")[0]


def load_map(path: str) -> dict:
    """Two columns: accession <tab> species-or-organism-text. Anything that is."""
    out = {}
    valid = set(SPECIES_COLUMNS) | {"RSV-unassigned"}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]", line, maxsplit=1)
            if len(parts) != 2:
                continue
            acc, text = parts[0].strip(), parts[1].strip()
            label = text if text in valid else classify_organism(text)
            if label:
                out[bare_accession(acc)] = label
    print(f"Loaded {len(out)} accession assignments from {path}")
    return out


def fetch_from_ncbi(accessions, email, cache_path, api_key=None,
                    batch=180) -> dict:
    """Resolve accessions via NCBI esummary. Cached, so a rerun costs nothing."""
    import urllib.parse
    import urllib.request

    cache = {}
    if cache_path and os.path.isfile(cache_path):
        try:
            cache = json.load(open(cache_path))
            print(f"Loaded {len(cache)} cached NCBI records from {cache_path}")
        except (OSError, ValueError):
            cache = {}

    todo = sorted({bare_accession(a) for a in accessions} - set(cache))
    if not todo:
        return cache
    print(f"Querying NCBI for {len(todo)} accessions...")

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        params = {"db": "nuccore", "id": ",".join(chunk),
                  "retmode": "json", "tool": "motif-species",
                  "email": email}
        if api_key:
            params["api_key"] = api_key
        try:
            with urllib.request.urlopen(
                    f"{url}?{urllib.parse.urlencode(params)}",
                    timeout=60) as resp:
                data = json.load(resp)
        except Exception as exc:
            print(f"  ! NCBI request failed: {exc}", file=sys.stderr)
            print("    Falling back to whatever is already resolved. Use "
                  "--map to supply the rest offline.", file=sys.stderr)
            break
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            acc = bare_accession(rec.get("accessionversion")
                                 or rec.get("caption") or "")
            title = " ".join(filter(None, [rec.get("organism", ""),
                                           rec.get("title", ""),
                                           rec.get("strain", "")]))
            if acc:
                cache[acc] = title.strip()
        time.sleep(0.4 if api_key else 1.0)   # NCBI rate limit

    if cache_path:
        try:
            with open(cache_path, "w") as fh:
                json.dump(cache, fh, indent=1, sort_keys=True)
            print(f"Cached {len(cache)} NCBI records to {cache_path}")
        except OSError as exc:
            print(f"  ! could not write cache: {exc}", file=sys.stderr)
    return cache


# Stockholm

def read_stockholm(path: str) -> dict:
    """Returns the raw lines plus the sequence names in order of appearance."""
    lines, names = [], []
    seen = set()
    with open(path, errors="replace") as fh:
        for line in fh:
            lines.append(line.rstrip("\n"))
            stripped = line.strip()
            if not stripped or stripped.startswith("#") \
                    or stripped.startswith("//"):
                continue
            name = stripped.split(None, 1)[0]
            if name not in seen:
                seen.add(name)
                names.append(name)
    return {"lines": lines, "names": names}


def write_annotated(sto, assignments, out_path, rename=False):
    """Write a copy carrying '#=GS <seq> DE <species>' for every sequence."""
    names = sto["names"]
    new_name = {n: f"{n}|{assignments.get(n, 'Unclassified')}" for n in names} \
        if rename else {n: n for n in names}

    gs_block = [f"#=GS {new_name[n]:<40s} DE {assignments.get(n, 'Unclassified')}"
                for n in names]

    # Stockholm convention: #=GF lines first, then #=GS, then the rows.
    out, inserted = [], False
    for line in sto["lines"]:
        stripped = line.strip()
        if (not inserted and stripped
                and not stripped.startswith("# STOCKHOLM")
                and not stripped.startswith("#=GF")):
            out.extend(gs_block)
            out.append("")
            inserted = True
        if rename and stripped and not stripped.startswith("//"):
            if stripped.startswith("#=GS ") or stripped.startswith("#=GR "):
                parts = line.split(None, 2)
                if len(parts) == 3 and parts[1] in new_name:
                    line = f"{parts[0]} {new_name[parts[1]]} {parts[2]}"
            elif not stripped.startswith("#"):
                parts = line.split(None, 1)
                if len(parts) == 2 and parts[0] in new_name:
                    line = f"{new_name[parts[0]]:<48s} {parts[1]}"
        out.append(line)
    if not inserted:
        out = gs_block + [""] + out

    with open(out_path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def collect_stockholms(root, series_filter):
    """{(series, motif): path} for every .stockholm under root."""
    found = {}
    for dirpath, _dirs, filenames in os.walk(root):
        for name in filenames:
            if not name.endswith(".stockholm"):
                continue
            series, num = parse_name(name)
            if num is None:
                continue
            if series_filter and series not in series_filter:
                continue
            key = (series, num)
            path = os.path.join(dirpath, name)
            if key in found:
                old_p, new_p = found[key], path
                keep = max((old_p, new_p), key=os.path.getmtime)
                drop = old_p if keep is new_p else new_p
                print(f"  ! {series} motif {num}: two .stockholm files\n"
                      f"        using   {os.path.relpath(keep, root)} "
                      f"({os.path.getsize(keep)} bytes)\n"
                      f"        ignored {os.path.relpath(drop, root)} "
                      f"({os.path.getsize(drop)} bytes)", file=sys.stderr)
                if os.path.getmtime(path) <= os.path.getmtime(found[key]):
                    continue
            found[key] = path
    return found


# The 29 covarying motifs, so the script runs with no CSV on hand.
# (gene, motif number, coordinates, alignment depth)
BUILTIN_MOTIFS = [
    ("N", 3, "477-575", 131),
    ("P", 14, "1577-1649", 81),
    ("P", 16, "1804-1916", 139),
    ("P", 17, "1985-2078", 131),
    ("M", 26, "2738-2797", 56),
    ("M", 27, "2809-2866", 52),
    ("F", 28, "3009-3092", 76),
    ("F", 30, "3370-3425", 33),
    ("F", 41, "4269-4371", 116),
    ("M2", 48, "4792-4811", 19),
    ("M2", 49, "4817-4830", 3),
    ("M2", 52, "4998-5012", 3),
    ("M2", 53, "5016-5074", 50),
    ("SH", 59, "5757-5823", 62),
    ("SH", 60, "5848-5935", 147),
    ("G", 65, "6149-6169", 5),
    ("G", 66, "6191-6254", 136),
    ("L", 69, "7093-7139", 116),
    ("L", 75, "7810-7844", 42),
    ("L", 79, "8314-8390", 95),
    ("L", 80, "8663-8749", 59),
    ("L", 82, "9054-9127", 90),
    ("L", 85, "9427-9480", 40),
    ("L", 91, "10470-10553", 121),
    ("L", 94, "10887-10947", 77),
    ("L", 96, "11261-11344", 106),
    ("L", 98, "11785-11863", 80),
    ("L", 101, "12128-12213", 102),
    ("L", 103, "12294-12392", 109),
]


def builtin_motif_list():
    return [{"motif": n, "gene": g, "coords": c, "expected_depth": d}
            for g, n, c, d in BUILTIN_MOTIFS]


def read_motif_list(path):
    """The motif set to process, from the covarying-motif CSV."""
    wanted = []
    with open(path, newline="", errors="replace") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        for row in csv.DictReader(fh, delimiter=delim):
            key = {k.lower().strip(): (v or "").strip()
                   for k, v in row.items() if k}
            num = None
            for k, v in key.items():
                if "motif" in k and v:
                    m = re.search(r"(\d+)", v)
                    if m:
                        num = int(m.group(1))
                        break
            if num is None:
                continue
            depth = key.get("alignment depth (n)") or ""
            wanted.append({
                "motif": num,
                "gene": key.get("gene", ""),
                "coords": key.get("coordinates", ""),
                "expected_depth": int(depth) if depth.isdigit() else None,
            })
    return wanted


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root",
                    default="/lustre/hdd/LAS/wmoss-lab/ddabiri/CaCo_folder/MNV")
    ap.add_argument("--motifs", default="builtin",
                    help="CSV/TSV of the motifs to process (needs Gene, "
                         "Motif no., Coordinates). Defaults to the built-in "
                         "list of the 29 covarying motifs; pass a file to "
                         "override, or 'all' for every motif found.")
    ap.add_argument("--series", action="append",
                    help="restrict to this ScanFold run, e.g. NC_hMNV_cRNA")
    ap.add_argument("--out", default="species_composition.csv")
    ap.add_argument("--annotated-dir", default="annotated_alignments")
    ap.add_argument("--audit", help="per-accession assignment and its source")
    ap.add_argument("--unresolved", default="unresolved_accessions.txt")
    ap.add_argument("--map", dest="map_file",
                    help="TSV of accession <tab> species or organism text")
    ap.add_argument("--fetch", action="store_true",
                    help="look unknown accessions up at NCBI")
    ap.add_argument("--email", help="contact address, required by NCBI")
    ap.add_argument("--api-key", help="NCBI API key (raises the rate limit)")
    ap.add_argument("--cache", default="ncbi_accessions.json")
    ap.add_argument("--rename", action="store_true",
                    help="also append the species to each sequence name")
    ap.add_argument("--no-annotate", action="store_true",
                    help="build the table only, write no alignments")
    args = ap.parse_args()

    if args.fetch and not args.email:
        sys.exit("--fetch requires --email (NCBI requires a contact address)")
    if not os.path.isdir(args.root):
        sys.exit(f"Not a directory: {args.root}")

    if args.motifs == "builtin":
        wanted = builtin_motif_list()
        print(f"{len(wanted)} motifs from the built-in covarying list")
    elif args.motifs == "all":
        wanted = None
        print("Processing every motif found on disk")
    elif not os.path.isfile(args.motifs):
        sys.exit(f"No such file: {args.motifs}. Pass a CSV, or use "
                 f"'builtin' (the 29 covarying motifs) or 'all'.")
    else:
        wanted = read_motif_list(args.motifs)
        print(f"{len(wanted)} motifs listed in {args.motifs}")

    series_filter = set(args.series) if args.series else None
    stockholms = collect_stockholms(args.root, series_filter)
    print(f"{len(stockholms)} .stockholm files found under {args.root}")

    if wanted is None:
        wanted = []
        for (series, num), path in sorted(stockholms.items(),
                                          key=lambda kv: kv[0][1]):
            m = COORD_RE.search(os.path.basename(path))
            wanted.append({"motif": num, "gene": "",
                           "coords": f"{m.group(1)}-{m.group(2)}" if m else "",
                           "expected_depth": None})

    # ---- gather sequence names -----------------------------------------
    per_motif, missing = {}, []
    for w in wanted:
        hits = [k for k in stockholms if k[1] == w["motif"]]
        if not hits:
            missing.append(w["motif"])
            continue
        if len(hits) > 1:
            print(f"  ! motif {w['motif']} matches {len(hits)} series; use "
                  f"--series", file=sys.stderr)
            continue
        path = stockholms[hits[0]]
        sto = read_stockholm(path)
        per_motif[w["motif"]] = {"w": w, "path": path, "sto": sto}
    if missing:
        print(f"  ! no .stockholm for motifs: "
              f"{', '.join(map(str, missing))}", file=sys.stderr)

    # ---- resolve every accession once ----------------------------------
    name_to_acc, all_accs = {}, set()
    for rec in per_motif.values():
        for n in rec["sto"]["names"]:
            m = ACCESSION_RE.search(n)
            acc = bare_accession(m.group(1)) if m else None
            name_to_acc[n] = acc
            if acc:
                all_accs.add(acc)
    print(f"{len(name_to_acc)} sequence names, "
          f"{len(all_accs)} distinct accessions")

    assign, source = {}, {}

    # 1. organism text already in the sequence name
    for n in name_to_acc:
        label = classify_organism(n)
        if label:
            assign[n], source[n] = label, "sequence name"

    # 2. user-supplied map
    acc_map = load_map(args.map_file) if args.map_file else {}

    # 3. built-ins
    acc_label = dict(BUILTIN_ACCESSIONS)
    acc_label.update(acc_map)          # the map overrides built-ins
    acc_src = {a: ("--map" if a in acc_map else "built-in table")
               for a in acc_label}

    for n, acc in name_to_acc.items():
        if n in assign or not acc or acc not in acc_label:
            continue
        assign[n], source[n] = acc_label[acc], acc_src[acc]

    # 4. NCBI for whatever is left
    unresolved_accs = {acc for n, acc in name_to_acc.items()
                       if n not in assign and acc}
    if unresolved_accs and args.fetch:
        titles = fetch_from_ncbi(unresolved_accs, args.email, args.cache,
                                 args.api_key)
        for n, acc in name_to_acc.items():
            if n in assign or not acc:
                continue
            label = classify_organism(titles.get(acc, ""))
            if label:
                assign[n], source[n] = label, "NCBI"
    elif unresolved_accs and os.path.isfile(args.cache):
        try:
            titles = json.load(open(args.cache))
            for n, acc in name_to_acc.items():
                if n in assign or not acc:
                    continue
                label = classify_organism(titles.get(acc, ""))
                if label:
                    assign[n], source[n] = label, "NCBI cache"
            print(f"Resolved from cache {args.cache}")
        except (OSError, ValueError):
            pass

    # When little or nothing resolved, the sequence-name format is usually
    # the reason. Show examples rather than leaving a table of zeros.
    n_unassigned = sum(1 for n in name_to_acc if n not in assign)
    if name_to_acc and n_unassigned > len(name_to_acc) // 2:
        examples = [n for n in sorted(name_to_acc) if n not in assign][:10]
        no_acc = sum(1 for n in name_to_acc
                     if n not in assign and not name_to_acc[n])
        print("\n" + "=" * 68, file=sys.stderr)
        print(f"{n_unassigned} of {len(name_to_acc)} sequences unassigned. "
              f"Example names:", file=sys.stderr)
        for e in examples:
            print(f"    {e}   -> accession parsed: "
                  f"{name_to_acc[e] or 'NONE'}", file=sys.stderr)
        if no_acc:
            print(f"\n  {no_acc} names contain no recognisable accession at "
                  f"all, so no lookup can help. Send me a few of these "
                  f"names and the pattern can be added.", file=sys.stderr)
        else:
            print(f"\n  Accessions parsed fine, so they just need "
                  f"resolving: rerun with --fetch --email <you@iastate.edu>",
                  file=sys.stderr)
        print("=" * 68 + "\n", file=sys.stderr)

    still = sorted({name_to_acc[n] or n for n in name_to_acc
                    if n not in assign})
    if still:
        with open(args.unresolved, "w") as fh:
            fh.write("\n".join(still) + "\n")
        print(f"\n! {len(still)} accessions could not be assigned to a "
              f"species. Written to {args.unresolved}.", file=sys.stderr)
        print(f"  They are counted under 'Unclassified', never guessed. "
              f"Resolve them with --fetch --email ..., or supply "
              f"--map accession_species.tsv.", file=sys.stderr)

    # ---- build the table ------------------------------------------------
    if not args.no_annotate:
        os.makedirs(args.annotated_dir, exist_ok=True)

    columns = list(SPECIES_COLUMNS)
    extra = sorted({v for v in assign.values() if v not in columns})
    rows, depth_warnings = [], []

    for w in wanted:
        rec = per_motif.get(w["motif"])
        if rec is None:
            continue
        names = rec["sto"]["names"]
        labels = [assign.get(n, "Unclassified") for n in names]
        counts = Counter(labels)

        row = {"Gene": w["gene"], "Motif No": w["motif"],
               "Coordinates": w["coords"]}
        for col in columns + extra:
            row[col] = counts.get(col, 0)
        unclassified = counts.get("Unclassified", 0)
        if unclassified:
            row["Unclassified"] = unclassified
        row["Total _Seq"] = len(names)
        rows.append(row)

        if w["expected_depth"] is not None and w["expected_depth"] != len(names):
            depth_warnings.append(
                f"motif {w['motif']}: {len(names)} sequences in the "
                f"alignment vs alignment depth {w['expected_depth']} in "
                f"{os.path.basename(args.motifs)}")

        if not args.no_annotate:
            base = os.path.basename(rec["path"])
            out_path = os.path.join(args.annotated_dir,
                                    f"{base}.species.stockholm")
            write_annotated(rec["sto"],
                            {n: assign.get(n, "Unclassified") for n in names},
                            out_path, args.rename)

    if not rows:
        sys.exit("No motifs could be processed.")

    cols = (["Gene", "Motif No", "Coordinates"] + columns + extra
            + (["Unclassified"] if any("Unclassified" in r for r in rows)
               else []) + ["Total _Seq"])
    for r in rows:
        for c in cols:
            r.setdefault(c, 0)

    delim = "\t" if args.out.lower().endswith((".tsv", ".txt")) else ","
    with open(args.out, "w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=cols, delimiter=delim,
                            lineterminator="\n")
        w_.writeheader()
        w_.writerows(rows)
    print(f"\nWrote {len(rows)} motifs to {args.out}")
    if not args.no_annotate:
        print(f"Wrote {len(rows)} annotated alignments to "
              f"{args.annotated_dir}/")

    if args.audit:
        with open(args.audit, "w", newline="") as fh:
            w_ = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w_.writerow(["Sequence name", "Accession", "Species", "Source"])
            for n in sorted(name_to_acc):
                w_.writerow([n, name_to_acc[n] or "NA",
                             assign.get(n, "Unclassified"),
                             source.get(n, "unresolved")])
        print(f"Wrote per-sequence assignments to {args.audit}")

    totals = Counter()
    for r in rows:
        for c in cols:
            if c not in ("Gene", "Motif No", "Coordinates"):
                totals[c] += r[c]
    print("\nSequences across all motifs:")
    for c in cols:
        if c in ("Gene", "Motif No", "Coordinates"):
            continue
        print(f"    {c:16s} {totals[c]}")
    by_source = Counter(source.get(n, "unresolved") for n in name_to_acc)
    print("\nAssignments by source:")
    for s, n in by_source.most_common():
        print(f"    {s:16s} {n}")

    if depth_warnings:
        print(f"\n! {len(depth_warnings)} depth mismatches:", file=sys.stderr)
        for d in depth_warnings:
            print(f"  ! {d}", file=sys.stderr)


if __name__ == "__main__":
    main()
