#!/usr/bin/env python3
"""Write motif_N.dbn and motif_N.sto for a list of motifs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import Counter

MOTIF_RE = re.compile(r"motif[_-](\d+)(?=[._\-]|$)", re.IGNORECASE)
COORD_RE = re.compile(r"coordinates[:_](\d+)[-_](\d+)", re.IGNORECASE)
ZSCORE_RE = re.compile(r"z[-_]?score[=:_](-?\d+(?:\.\d+)?)", re.IGNORECASE)
ACCESSION_RE = re.compile(r"\b([A-Z]{1,2}_?\d{5,8}(?:\.\d+)?)\b")

GAPS = set("-.~_")
PAIRED_OPEN = "([{<"
PAIRED_CLOSE = ")]}>"
PAIRED_CH = set(PAIRED_OPEN + PAIRED_CLOSE + "AaBbCcDd")
LOOP_CH = set(".,:_-~")
ALN_EXT = (".stockholm", ".sto", ".stk", ".fa", ".fasta", ".fna", ".afa",
           ".aln")

# Suffixes that mark a derived file rather than the primary alignment.
DERIVED = (".bak", ".ps", ".svg", ".pdf", ".png", ".eps")


# Default motif set. Replace, extend, or override with --motifs FILE;
# --all-motifs exports everything found instead.
DEFAULT_MOTIFS = [
    3, 14, 16, 17, 26, 27, 28, 30, 41, 48, 53, 59, 60, 66, 69, 75, 79, 80,
    82, 85, 91, 94, 96, 98, 101, 103,
]


def parse_name(basename):
    m = MOTIF_RE.search(basename)
    if not m:
        return None, None
    series = basename[:m.start()].rstrip("_-.") or "(unnamed)"
    return series, int(m.group(1))


# Readers

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


def read_fasta(path):
    names, seqs, cur = [], [], []
    for line in open(path, encoding="utf-8", errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            names.append(line[1:].strip())
            cur = []
        elif line:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return names, [s.upper().replace("T", "U") for s in seqs], {}


def read_stockholm(path):
    """Interleaved blocks are merged: a repeated name continues one row."""
    names, seqs, gc, index = [], [], {}, {}
    for line in open(path, encoding="utf-8", errors="ignore"):
        s = line.rstrip("\n").strip()
        if not s or s == "//":
            continue
        if s.startswith("#=GC"):
            p = s.split(None, 2)
            if len(p) == 3:
                gc[p[1]] = gc.get(p[1], "") + p[2].strip()
            continue
        if s.startswith("#"):
            continue
        p = s.rsplit(None, 1)
        if len(p) != 2:
            continue
        nm, sq = p[0].strip(), p[1].strip()
        if nm in index:
            seqs[index[nm]] += sq
        else:
            index[nm] = len(names)
            names.append(nm)
            seqs.append(sq)
    return names, [s.upper().replace("T", "U") for s in seqs], gc


def read_alignment(path):
    return read_fasta(path) if sniff(path) == "fasta" else read_stockholm(path)


def read_dbn_file(path):
    """Last dot-bracket line in a .dbn, ignoring header and sequence."""
    for line in reversed(open(path, encoding="utf-8",
                              errors="ignore").read().splitlines()):
        tok = line.strip().split()[0] if line.strip() else ""
        if not tok or tok.startswith(">"):
            continue
        if set(tok) <= (PAIRED_CH | LOOP_CH) and any(c in "()[]{}<>"
                                                     for c in tok):
            return tok
    return None


def read_pairs_from_cov(path, length):
    """Rebuild a dot-bracket from the pair list in an R-scape .cov/.power."""
    pairs = []
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("#") or not line.strip():
            continue
        parts = line.replace("*", " ").split()
        if len(parts) < 2:
            continue
        try:
            i, j = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        if 1 <= i < j <= length:
            pairs.append((i, j))
    if not pairs:
        return None, 0
    pairs.sort()
    used, kept, dropped = set(), [], 0
    for i, j in pairs:
        if i in used or j in used:
            dropped += 1
            continue
        if any(a < i < b < j or i < a < j < b for a, b in kept):
            dropped += 1
            continue
        kept.append((i, j))
        used.update((i, j))
    out = ["."] * length
    for i, j in kept:
        out[i - 1], out[j - 1] = "(", ")"
    return "".join(out), dropped


# Structure helpers

def normalise_structure(ss, keep_pseudoknots=False):
    """Collapse WUSS / VARNA notation to plain dot-bracket."""
    out = []
    for ch in ss:
        if ch in "(<":
            out.append("(")
        elif ch in ")>":
            out.append(")")
        elif ch in "[{":
            out.append(ch if keep_pseudoknots else ".")
        elif ch in "]}":
            out.append(ch if keep_pseudoknots else ".")
        elif ch in "ABCD":
            out.append("[" if keep_pseudoknots else ".")
        elif ch in "abcd":
            out.append("]" if keep_pseudoknots else ".")
        else:
            out.append(".")
    return "".join(out)


def rebalance(ss):
    """Remove brackets left without a partner."""
    chars = list(ss)
    removed = 0
    for o, c in (("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")):
        stack, orphans = [], []
        for i, ch in enumerate(chars):
            if ch == o:
                stack.append(i)
            elif ch == c:
                if stack:
                    stack.pop()
                else:
                    orphans.append(i)
        for i in stack + orphans:
            chars[i] = "."
            removed += 1
    return "".join(chars), removed


def check_balanced(ss):
    """Return a list of problems; empty means the structure is well formed."""
    problems = []
    for o, c in zip("([{<", ")]}>"):
        stack = 0
        for ch in ss:
            if ch == o:
                stack += 1
            elif ch == c:
                stack -= 1
                if stack < 0:
                    problems.append(f"unmatched '{c}'")
                    break
        if stack > 0:
            problems.append(f"{stack} unmatched '{o}'")
    return problems


def degap_structure(ss, row):
    """Structure restricted to the non-gap columns of one sequence."""
    return "".join(ss[i] for i, ch in enumerate(row) if ch not in GAPS)


def structure_for_columns(dbn, seqs, names, gc, ref_name=None):
    """Place an ungapped .dbn onto gapped alignment columns."""
    L = len(seqs[0])
    if len(dbn) == L:
        return dbn, "structure already spans every column"

    def place(row, label):
        cols = [i for i, c in enumerate(row) if c not in GAPS]
        if len(cols) != len(dbn):
            return None
        out = ["."] * L
        for k, i in enumerate(cols):
            out[i] = dbn[k]
        return "".join(out), label

    if ref_name:
        for nm, sq in zip(names, seqs):
            if ref_name.lower() in nm.lower():
                r = place(sq, f"mapped via {nm!r}")
                if r:
                    return r
    rf = gc.get("RF")
    if rf and len(rf) == L:
        r = place(rf, "mapped via #=GC RF")
        if r:
            return r
    for idx, sq in enumerate(seqs):
        r = place(sq, f"mapped via {names[idx]!r}")
        if r:
            return r
    return None, None


# Collecting inputs

def collect(root, series_filter, recursive=True):
    """{(series, motif): {'aln': path, 'dbn': path, 'cov': path}}"""
    found = {}
    walker = os.walk(root) if recursive else [
        (root, [], [f for f in os.listdir(root)
                    if os.path.isfile(os.path.join(root, f))])]
    for dirpath, _dirs, filenames in walker:
        for name in filenames:
            low = name.lower()
            if low.endswith(DERIVED):
                continue
            series, num = parse_name(name)
            if num is None:
                continue
            if series_filter and series not in series_filter:
                continue
            path = os.path.join(dirpath, name)
            rec = found.setdefault((series, num), {})

            if low.endswith(".dbn"):
                slot = "dbn"
            elif low.endswith(ALN_EXT) and ".r2r." not in low:
                # Prefer the plain alignment over a CaCoFold rewrite unless
                # the user asks for the latter.
                slot = "aln_cacofold" if "cacofold" in low else "aln"
            elif low.endswith(".power") or low.endswith(".cov"):
                slot = "cov_cacofold" if "cacofold" in low else "cov"
            else:
                continue

            if slot in rec:
                # keep the newer file, and say so
                if os.path.getmtime(path) <= os.path.getmtime(rec[slot]):
                    continue
                rec.setdefault("_dupes", []).append(
                    (slot, rec[slot], path))
            rec[slot] = path
    return found


def read_motif_list(path):
    nums = []
    with open(path, newline="", errors="replace") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        rows = list(csv.DictReader(fh, delimiter=delim))
        if rows and rows[0]:
            for row in rows:
                low = {(k or "").strip().lower(): (v or "").strip()
                       for k, v in row.items() if k}
                for key, val in low.items():
                    if "motif" in key and val:
                        m = re.search(r"(\d+)", val)
                        if m:
                            nums.append(int(m.group(1)))
                            break
    if not nums:          # plain list of numbers, one per line
        for line in open(path, errors="replace"):
            m = re.search(r"(\d+)", line)
            if m:
                nums.append(int(m.group(1)))
    return sorted(set(nums))


def load_name_map(path):
    """accession -> label, from a two-column TSV/CSV."""
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


def display_name(raw, name_map, style, keep_version=True):
    """Sequence header for the output file."""
    m = ACCESSION_RE.search(raw)
    acc_full = m.group(1) if m else raw.split("/")[0].split()[0]
    acc_bare = acc_full.split(".")[0]
    acc = acc_full if keep_version else acc_bare
    if style == "raw":
        return raw
    label = name_map.get(acc_bare)
    if style == "accession" or not label:
        return acc
    if style == "species":
        return label
    return f"{acc}_{label}"          # 'both'


# Writers

def wrap(seq, width):
    if width <= 0:
        return [seq]
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def write_dbn(path, structure, seq=None, header=None, width=0):
    with open(path, "w") as fh:
        if header:
            fh.write(f">{header}\n")
        if seq is not None:
            for chunk in wrap(seq, width):
                fh.write(chunk + "\n")
        for chunk in wrap(structure, width):
            fh.write(chunk + "\n")


def write_fasta_sto(path, headers, seqs, width=0):
    with open(path, "w") as fh:
        for h, s in zip(headers, seqs):
            fh.write(f">{h}\n")
            for chunk in wrap(s, width):
                fh.write(chunk + "\n")


def write_stockholm(path, headers, seqs, structure=None):
    pad = max((len(h) for h in headers), default=10)
    pad = max(pad, len("#=GC SS_cons"))
    with open(path, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n\n")
        for h, s in zip(headers, seqs):
            fh.write(f"{h:<{pad}} {s}\n")
        if structure:
            fh.write(f"{'#=GC SS_cons':<{pad}} {structure}\n")
        fh.write("//\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".",
                    help="folder holding the alignments and structures "
                         "(default: the current directory)")
    ap.add_argument("--out", default="exported_motifs",
                    help="output directory (default exported_motifs)")
    ap.add_argument("--series", action="append",
                    help="restrict to one ScanFold run, e.g. NC_hMNV_cRNA")
    ap.add_argument("--motifs",
                    help="CSV/TSV or plain list of motif numbers to export, "
                         "overriding the built-in default set")
    ap.add_argument("--all-motifs", action="store_true",
                    help="export every motif found, ignoring the built-in "
                         "default set")
    ap.add_argument("--prefix", default="motif_",
                    help="output filename prefix (default 'motif_', giving "
                         "motif_13.dbn and motif_13.sto)")

    ap.add_argument("--sto-format", choices=["fasta", "stockholm", "both"],
                    default="fasta",
                    help="fasta (default) writes the >accession layout into "
                         ".sto; stockholm writes real Stockholm; both writes "
                         ".sto and .stockholm")
    ap.add_argument("--use-cacofold", action="store_true",
                    help="prefer the .cacofold alignment and structure")
    ap.add_argument("--ss-from-cov", action="store_true",
                    help="if no .dbn and no SS_cons, rebuild the structure "
                         "from the R-scape pair list")
    ap.add_argument("--keep-pseudoknots", action="store_true",
                    help="keep crossing pairs as [] brackets; by default "
                         "they are written as unpaired, since plain "
                         "dot-bracket cannot represent them")
    ap.add_argument("--ref",
                    help="substring of the sequence the .dbn describes, used "
                         "when mapping it onto alignment columns")

    ap.add_argument("--degap-to",
                    help="write the structure degapped to this sequence "
                         "(substring match) instead of to alignment columns")
    ap.add_argument("--dbn-with-sequence", action="store_true",
                    help="include the reference sequence line in the .dbn")
    ap.add_argument("--dbn-header", action="store_true",
                    help="write a '>motif_N' header line into the .dbn")

    ap.add_argument("--names", choices=["accession", "species", "both", "raw"],
                    default="accession",
                    help="sequence header style (default accession)")
    ap.add_argument("--name-map",
                    help="two-column accession->label file for --names")
    ap.add_argument("--no-version", action="store_true",
                    help="strip the .1 version suffix from accessions")
    ap.add_argument("--wrap", type=int, default=0,
                    help="wrap sequences at N characters (0 = one line)")
    ap.add_argument("--min-seqs", type=int, default=1,
                    help="skip motifs with fewer sequences than this")
    ap.add_argument("--manifest", default="export_manifest.tsv",
                    help="per-motif record of what was written")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be written, create nothing")
    args = ap.parse_args()

    if not os.path.isdir(args.root):
        sys.exit(f"Not a directory: {args.root}")

    series_filter = set(args.series) if args.series else None
    found = collect(args.root, series_filter)
    if not found:
        sys.exit(f"No motif files found under {args.root}")

    counts = Counter(s for s, _n in found)
    print(f"Found {len(found)} motifs under {args.root}")
    for s, n in sorted(counts.items()):
        print(f"    series {s!r}: {n} motifs")
    if len(counts) > 1 and not series_filter:
        print("\n!! More than one run is present and each numbers its motifs "
              "independently, so motif_4 of one is not motif_4 of the other. "
              "Use --series to pick one.\n", file=sys.stderr)

    wanted, where = None, ""
    if args.all_motifs:
        print("Exporting every motif found (--all-motifs)")
    elif args.motifs:
        if not os.path.isfile(args.motifs):
            sys.exit(f"No such file: {args.motifs}")
        wanted = set(read_motif_list(args.motifs))
        where = f"from {args.motifs}"
    else:
        wanted = set(DEFAULT_MOTIFS)
        where = "from the built-in default set (--all-motifs for every motif)"
    if wanted is not None:
        print(f"Exporting {len(wanted)} motifs {where}: "
              f"{', '.join(str(n) for n in sorted(wanted))}")
        absent = sorted(wanted - {n for _s, n in found})
        if absent:
            print(f"  ! requested but not found on disk: "
                  f"{', '.join(map(str, absent))}", file=sys.stderr)

    name_map = load_name_map(args.name_map) if args.name_map else {}
    if args.name_map:
        print(f"Loaded {len(name_map)} name mappings")
    if args.names in ("species", "both") and not name_map:
        print("  ! --names needs --name-map to supply labels; falling back "
              "to accessions", file=sys.stderr)

    if not args.dry_run:
        os.makedirs(args.out, exist_ok=True)

    manifest, skipped = [], []
    for (series, num) in sorted(found, key=lambda k: (k[0], k[1])):
        rec = found[(series, num)]
        if wanted is not None and num not in wanted:
            continue

        aln_key = ("aln_cacofold" if args.use_cacofold and
                   "aln_cacofold" in rec else "aln")
        aln_path = rec.get(aln_key) or rec.get("aln") or rec.get("aln_cacofold")
        if not aln_path:
            skipped.append((num, "no alignment file"))
            continue

        try:
            names, seqs, gc = read_alignment(aln_path)
        except OSError as exc:
            skipped.append((num, f"cannot read alignment: {exc}"))
            continue
        if not seqs:
            skipped.append((num, "alignment has no sequences"))
            continue
        if len(seqs) < args.min_seqs:
            skipped.append((num, f"only {len(seqs)} sequences"))
            continue

        L = len(seqs[0])
        ragged = sorted({len(s) for s in seqs})
        if len(ragged) > 1:
            skipped.append((num, f"rows differ in length {ragged[:5]}"))
            continue

        # ---- structure -------------------------------------------------
        source, structure = None, None
        dbn_path = rec.get("dbn")
        if dbn_path:
            raw = read_dbn_file(dbn_path)
            if raw:
                placed, how = structure_for_columns(raw, seqs, names, gc,
                                                    args.ref)
                if placed:
                    structure, source = placed, f"{os.path.basename(dbn_path)} ({how})"
                else:
                    ung = sorted({sum(1 for c in s if c not in GAPS)
                                  for s in seqs})
                    skipped.append(
                        (num, f".dbn is {len(raw)} nt but the alignment has "
                              f"{L} columns (ungapped {ung[:4]}); cannot map"))
                    continue
        if structure is None:
            ss = gc.get("SS_cons")
            if ss and len(ss) == L:
                structure, source = ss, "#=GC SS_cons"
        if structure is None and args.ss_from_cov:
            cov = (rec.get("cov_cacofold") if args.use_cacofold
                   else rec.get("cov")) or rec.get("cov") \
                or rec.get("cov_cacofold")
            if cov:
                rebuilt, dropped = read_pairs_from_cov(cov, L)
                if rebuilt:
                    structure = rebuilt
                    source = f"{os.path.basename(cov)} pair list"
                    if dropped:
                        source += f", {dropped} crossing pair(s) dropped"
        if structure is None:
            skipped.append((num, "no .dbn, no SS_cons, no pair list"))
            continue

        structure = normalise_structure(structure,
                                        args.keep_pseudoknots)
        problems = check_balanced(structure)
        if problems:
            print(f"  ! motif {num}: structure is not balanced "
                  f"({'; '.join(problems)})", file=sys.stderr)

        # ---- optional degapping to one sequence -------------------------
        out_seqs, out_structure = list(seqs), structure
        degapped_to = ""
        if args.degap_to:
            hit = next((i for i, nm in enumerate(names)
                        if args.degap_to.lower() in nm.lower()), None)
            if hit is None:
                print(f"  ! motif {num}: --degap-to {args.degap_to!r} matches "
                      f"no sequence; keeping all columns", file=sys.stderr)
            else:
                keep = [i for i, ch in enumerate(seqs[hit]) if ch not in GAPS]
                out_seqs = ["".join(s[i] for i in keep) for s in seqs]
                out_structure = "".join(structure[i] for i in keep)
                degapped_to = names[hit]

        headers = [display_name(n, name_map, args.names,
                                not args.no_version) for n in names]
        dupes = [h for h, c in Counter(headers).items() if c > 1]
        if dupes:
            print(f"  ! motif {num}: duplicate headers after renaming "
                  f"({', '.join(dupes[:3])}); appending an index",
                  file=sys.stderr)
            seen = Counter()
            fixed = []
            for h in headers:
                seen[h] += 1
                fixed.append(h if seen[h] == 1 else f"{h}_{seen[h]}")
            headers = fixed

        stem = f"{args.prefix}{num}"
        dbn_out = os.path.join(args.out, stem + ".dbn")
        sto_out = os.path.join(args.out, stem + ".sto")
        stk_out = os.path.join(args.out, stem + ".stockholm")

        # The .dbn carries the ungapped structure for the reference row
        # unless the alignment was already degapped.
        if args.degap_to and degapped_to:
            dbn_structure = out_structure
            dbn_seq = out_seqs[names.index(degapped_to)]
        else:
            # The .dbn must describe one ungapped sequence. Prefer the row
            # the user named, else the row with the fewest gaps -- degapping
            # to a gappy row would delete half of real base pairs.
            ref_i = next((i for i, nm in enumerate(names)
                          if args.ref and args.ref.lower() in nm.lower()),
                         None)
            if ref_i is None:
                ref_i = min(range(len(seqs)),
                            key=lambda i: sum(1 for c in seqs[i]
                                              if c in GAPS))
            dbn_structure = degap_structure(structure, seqs[ref_i])
            dbn_seq = "".join(c for c in seqs[ref_i] if c not in GAPS)

        dbn_structure, orphaned = rebalance(dbn_structure)
        if orphaned:
            print(f"  ! motif {num}: {orphaned} base pair(s) lost a partner "
                  f"when the structure was degapped to "
                  f"{names[ref_i] if not degapped_to else degapped_to!r}; "
                  f"those positions are now unpaired", file=sys.stderr)
        out_structure, _ = rebalance(out_structure)

        written = []
        if not args.dry_run:
            write_dbn(dbn_out, dbn_structure,
                      dbn_seq if args.dbn_with_sequence else None,
                      stem if args.dbn_header else None, args.wrap)
            written.append(os.path.basename(dbn_out))
            if args.sto_format in ("fasta", "both"):
                write_fasta_sto(sto_out, headers, out_seqs, args.wrap)
                written.append(os.path.basename(sto_out))
            if args.sto_format in ("stockholm", "both"):
                target = stk_out if args.sto_format == "both" else sto_out
                write_stockholm(target, headers, out_seqs, out_structure)
                written.append(os.path.basename(target))

        npairs = sum(1 for c in dbn_structure if c in PAIRED_OPEN)
        manifest.append({
            "motif": num, "series": series,
            "n_seq": len(out_seqs),
            "aln_columns": L,
            "dbn_length": len(dbn_structure),
            "base_pairs": npairs,
            "structure_source": source,
            "degapped_to": degapped_to or "alignment columns",
            "alignment_file": os.path.relpath(aln_path, args.root),
            "files": ", ".join(written) if written else "(dry run)",
        })

    if not manifest:
        if skipped:
            print(f"\n{len(skipped)} motif(s) skipped:", file=sys.stderr)
            for num, why in skipped:
                print(f"  ! motif {num}: {why}", file=sys.stderr)
        sys.exit("\nNothing was exported. Every motif lacked a usable "
                 "structure or alignment; the reasons are above. "
                 "--ss-from-cov can rebuild a structure from the R-scape "
                 "pair list when there is no .dbn and no SS_cons.")

    verb = "Would write" if args.dry_run else "Wrote"
    n_files = sum(len(m["files"].split(", ")) for m in manifest
                  if m["files"] != "(dry run)")
    print(f"\n{verb} {len(manifest)} motif(s)"
          + (f", {n_files} files, to {args.out}/" if not args.dry_run else ""))

    if not args.dry_run and args.manifest:
        mpath = os.path.join(args.out, args.manifest)
        with open(mpath, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(manifest[0]),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(manifest)
        print(f"Wrote {mpath}")

    lengths = Counter(m["dbn_length"] for m in manifest)
    nopairs = [m["motif"] for m in manifest if m["base_pairs"] == 0]
    print(f"\nStructure lengths range "
          f"{min(lengths)}-{max(lengths)} nt; "
          f"{sum(m['base_pairs'] for m in manifest)} base pairs in total")
    if nopairs:
        print(f"  ! {len(nopairs)} motif(s) have a structure with no base "
              f"pairs at all: {', '.join(map(str, nopairs))}",
              file=sys.stderr)
    if skipped:
        print(f"\n{len(skipped)} motif(s) skipped:", file=sys.stderr)
        for num, why in skipped:
            print(f"  ! motif {num}: {why}", file=sys.stderr)

    if args.sto_format == "fasta":
        print("\nNote: the .sto files are aligned FASTA, as requested. "
              "Infernal, R-scape and Easel will not read them; use "
              "--sto-format stockholm or both if you need that.")


if __name__ == "__main__":
    main()
