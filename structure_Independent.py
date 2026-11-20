#!/usr/bin/env python3
"""Whole-genome structure-agnostic covariation test. Superseded by pergene_agnostic.py."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import Counter

GAPS = set("-.~_")
PAIRED_CH = set("()[]{}<>AaBbCcDd")
LOOP_CH = set(".,:_-~")
MOTIF_RE = re.compile(r"motif[_-](\d+)(?=[._\-]|$)", re.IGNORECASE)
COORD_RE = re.compile(r"(\d+)\s*[-–]\s*(\d+)")

MAFFT_PRESETS = {
    "fast": ["--retree", "2", "--maxiterate", "0"],       # FFT-NS-2
    "auto": ["--auto"],
    "accurate": ["--localpair", "--maxiterate", "1000"],  # L-INS-i
}


# IO

def read_fasta(path):
    names, seqs, cur = [], [], None
    for line in open(path, encoding="utf-8", errors="ignore"):
        line = line.rstrip()
        if line.startswith(">"):
            names.append(line[1:].strip())
            seqs.append([])
            cur = len(seqs) - 1
        elif line and cur is not None:
            seqs[cur].append(line.strip())
    return names, ["".join(s).upper().replace("T", "U") for s in seqs]


def write_fasta(path, names, seqs):
    with open(path, "w") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def sto_names(names):
    """Stockholm sequence names cannot contain whitespace, and an NCBI FASTA."""
    out, seen = [], {}
    for n in names:
        base = n.split()[0] if n.split() else "seq"
        base = re.sub(r"[^A-Za-z0-9_.:/|-]", "_", base)
        seen[base] = seen.get(base, 0) + 1
        out.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return out


def write_sto(path, names, seqs, ss=None):
    names = sto_names(names)
    pad = max([len(n) for n in names] + [12]) + 2
    with open(path, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n\n")
        for n, s in zip(names, seqs):
            fh.write(f"{n:<{pad}}{s}\n")
        if ss:
            fh.write(f"{'#=GC SS_cons':<{pad}}{ss}\n")
        fh.write("//\n")


def read_dbn(path):
    for line in reversed(open(path, encoding="utf-8",
                              errors="ignore").read().splitlines()):
        tok = line.strip().split()[0] if line.strip() else ""
        if not tok or tok.startswith(">"):
            continue
        if set(tok) <= (PAIRED_CH | LOOP_CH) and any(c in "()[]{}<>"
                                                     for c in tok):
            return tok
    return None


def read_table(path):
    with open(path, newline="", errors="replace") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        d = "\t" if sample.count("\t") > sample.count(",") else ","
        return list(csv.DictReader(fh, delimiter=d))


def read_motif_coords(path):
    """{motif: (start, end)} from any table carrying motif and coordinates."""
    out = {}
    for row in read_table(path):
        low = {(k or "").strip().lower(): (v or "").strip()
               for k, v in row.items() if k}
        num = None
        for k, v in low.items():
            if "motif" in k and v:
                m = re.search(r"(\d+)", v)
                if m:
                    num = int(m.group(1))
                    break
        if num is None:
            continue
        s = low.get("start")
        e = low.get("end")
        if not (s and e):
            for k, v in low.items():
                if "coord" in k and v:
                    m = COORD_RE.search(v)
                    if m:
                        s, e = m.group(1), m.group(2)
                        break
        if s and e:
            out[num] = (int(s), int(e))
    return out


# Structure

def normalise(ss):
    return "".join("(" if c in "(<" else ")" if c in ")>" else "." for c in ss)


def pairs_of(ss):
    st, pr = [], set()
    for i, c in enumerate(ss):
        if c == "(":
            st.append(i)
        elif c == ")" and st:
            pr.add((st.pop(), i))
    return pr


def rebalance(ss):
    ch = list(ss)
    st, orph = [], []
    for i, c in enumerate(ch):
        if c == "(":
            st.append(i)
        elif c == ")":
            st.pop() if st else orph.append(i)
    for i in st + orph:
        ch[i] = "."
    return "".join(ch)


# Coordinate mapping

def coord_to_column(ref_row):
    """Map 1-based reference genome positions to alignment column indices."""
    return [i for i, c in enumerate(ref_row) if c not in GAPS]


def find_reference(names, token):
    """Index of the row whose name contains the reference token."""
    t = token.lower()
    exact = [i for i, n in enumerate(names) if n.lower().startswith(t)]
    if exact:
        return exact[0]
    loose = [i for i, n in enumerate(names) if t in n.lower()]
    return loose[0] if loose else None


def slice_columns(seqs, cols):
    return ["".join(s[c] for c in cols) for s in seqs]


def drop_empty_rows(names, seqs, min_residues=1):
    """Remove rows that are all gaps in the extracted region."""
    keep = [i for i, s in enumerate(seqs)
            if sum(1 for c in s if c not in GAPS) >= min_residues]
    return ([names[i] for i in keep], [seqs[i] for i in keep],
            len(names) - len(keep))


def drop_empty_columns(seqs, ss):
    """Remove columns that are gaps in every row, keeping ss in register."""
    L = len(seqs[0])
    keep = [i for i in range(L)
            if any(s[i] not in GAPS for s in seqs)]
    if len(keep) == L:
        return seqs, ss, 0
    return ([("".join(s[i] for i in keep)) for s in seqs],
            rebalance("".join(ss[i] for i in keep)),
            L - len(keep))


# Tools

def run(cmd, dry, log=None, stdout_to=None):
    if dry:
        print(f"    $ {' '.join(cmd)}"
              + (f" > {stdout_to}" if stdout_to else ""))
        return 0
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except FileNotFoundError:
        print(f"    ! not found: {cmd[0]}", file=sys.stderr)
        return 127
    except subprocess.TimeoutExpired:
        print(f"    ! timed out: {cmd[0]}", file=sys.stderr)
        return 124
    if log:
        with open(log, "w") as fh:
            fh.write(" ".join(cmd) + "\n" + p.stdout + p.stderr)
    if stdout_to and p.returncode == 0:
        with open(stdout_to, "w") as fh:
            fh.write(p.stdout)
    if p.returncode != 0:
        print(f"    ! exit {p.returncode}: {' '.join(cmd[:3])} ...",
              file=sys.stderr)
        for t in (p.stderr or p.stdout).strip().splitlines()[-3:]:
            print(f"      {t}", file=sys.stderr)
    return p.returncode


OBS = re.compile(r"BPAIRS\s+observed\s+to\s+covary\s*[:=]?\s*(\d+)", re.I)
EXP = re.compile(r"BPAIRS\s+expected\s+to\s+covary\s*[:=]?\s*(\d+(?:\.\d+)?)",
                 re.I)


def read_rscape(outdir):
    if not os.path.isdir(outdir):
        return None, None
    files = sorted(f for f in os.listdir(outdir)
                   if f.endswith((".cov", ".power")))
    files.sort(key=lambda f: ("cacofold" in f.lower(), f))
    for f in files:
        try:
            text = open(os.path.join(outdir, f), errors="replace").read()
        except OSError:
            continue
        m = OBS.search(text)
        if m:
            e = EXP.search(text)
            return int(m.group(1)), (float(e.group(1)) if e else None)
    return None, None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True,
                    help="FASTA of related genomes, including the reference")
    ap.add_argument("--reference", required=True,
                    help="accession or name substring of the genome whose "
                         "coordinates the motif table uses, e.g. NC_039199")
    ap.add_argument("--motifs", required=True,
                    help="table with motif numbers and coordinates")
    ap.add_argument("--structures",
                    help="folder of per-motif .dbn files. Without it the "
                         "regions are extracted but no structure is tested.")
    ap.add_argument("--out", default="agnostic")
    ap.add_argument("--alignment",
                    help="use this existing genome alignment instead of "
                         "running MAFFT")
    ap.add_argument("--mafft-preset", choices=list(MAFFT_PRESETS),
                    default="auto")
    ap.add_argument("--mafft-extra", default="--quiet --thread -1")
    ap.add_argument("--max-genomes", type=int,
                    help="subsample to this many genomes (the reference is "
                         "always kept). Use for a trial run.")
    ap.add_argument("--flank", type=int, default=0,
                    help="extra columns either side of each motif")
    ap.add_argument("--min-residues", type=int, default=10,
                    help="drop a sequence from a motif region if it has "
                         "fewer than this many non-gap residues there")
    ap.add_argument("--rscape-bin", default="R-scape",
                    help="path to the R-scape executable if it is not "
                         "simply 'R-scape' on PATH")
    ap.add_argument("--mafft-bin", default="mafft",
                    help="path to the MAFFT executable")
    ap.add_argument("--rscape-args", default="-s --nofigures",
                    help="R-scape flags. Use the SAME settings as your "
                         "original analysis, otherwise a difference in the "
                         "covariation counts could come from the settings "
                         "rather than from the alignment. Drop --nofigures "
                         "to have R-scape generate R2R figures.")
    ap.add_argument("--compare",
                    help="earlier results table to place beside the new "
                         "counts (needs Motif and a covariation column)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for p in (args.genomes, args.motifs):
        if not os.path.isfile(p):
            sys.exit(f"No such file: {p}")
    coords = read_motif_coords(args.motifs)
    if not coords:
        sys.exit(f"No motif coordinates could be read from {args.motifs}")
    print(f"{len(coords)} motif coordinate(s) from {args.motifs}")

    if not args.dry_run:
        need = ([args.rscape_bin]
                + ([] if args.alignment else [args.mafft_bin]))
        miss = [t for t in need
                if not (shutil.which(t) or os.path.isfile(t))]
        if miss:
            sys.exit(f"Not on PATH: {', '.join(miss)}")

    os.makedirs(args.out, exist_ok=True)

    # ---- 1-2. genome alignment ------------------------------------------
    if args.alignment:
        aln_path = args.alignment
        print(f"Using existing alignment {aln_path}")
    else:
        names, seqs = read_fasta(args.genomes)
        print(f"{len(names)} genome(s), lengths "
              f"{min(len(s) for s in seqs)}-{max(len(s) for s in seqs)} nt")
        ref_i = find_reference(names, args.reference)
        if ref_i is None:
            sys.exit(f"Reference {args.reference!r} is not in {args.genomes}. "
                     f"First few names: {', '.join(names[:5])}")
        print(f"Reference: {names[ref_i]} ({len(seqs[ref_i])} nt)")

        if args.max_genomes and len(names) > args.max_genomes:
            keep = [ref_i] + [i for i in range(len(names))
                              if i != ref_i][:args.max_genomes - 1]
            keep.sort()
            names = [names[i] for i in keep]
            seqs = [seqs[i] for i in keep]
            print(f"  subsampled to {len(names)} genomes "
                  f"(reference retained)")

        infa = os.path.join(args.out, "genomes_input.fa")
        aln_path = os.path.join(args.out, "genomes_aligned.fa")
        if not args.dry_run:
            write_fasta(infa, names, seqs)
        cmd = ([args.mafft_bin] + MAFFT_PRESETS[args.mafft_preset]
               + args.mafft_extra.split() + [infa])
        print(f"\nAligning with MAFFT ({args.mafft_preset}). This is the "
              f"slow step.")
        if run(cmd, args.dry_run, os.path.join(args.out, "mafft.log"),
               stdout_to=aln_path) != 0 and not args.dry_run:
            sys.exit("MAFFT failed; see mafft.log")

    if args.dry_run:
        print("\n(dry run: the remaining steps need the alignment)")
        return

    # ---- 3. coordinate map ----------------------------------------------
    anames, aseqs = read_fasta(aln_path)
    if not aseqs:
        sys.exit(f"No sequences in {aln_path}")
    if len({len(s) for s in aseqs}) != 1:
        sys.exit("Alignment rows differ in length")
    ncol = len(aseqs[0])
    ref_i = find_reference(anames, args.reference)
    if ref_i is None:
        sys.exit(f"Reference {args.reference!r} not in the alignment")
    colmap = coord_to_column(aseqs[ref_i])
    reflen = len(colmap)
    print(f"\nAlignment: {len(anames)} sequences, {ncol} columns")
    print(f"Reference {anames[ref_i]} spans {reflen} nt "
          f"({ncol - reflen} insertion columns in other genomes)")

    # ---- 4-6. per motif --------------------------------------------------
    rows = []
    for num in sorted(coords):
        start, end = coords[num]
        if end > reflen:
            print(f"motif {num}: coordinates {start}-{end} exceed the "
                  f"reference length {reflen}, skipped", file=sys.stderr)
            continue
        lo = max(0, start - 1 - args.flank)
        hi = min(reflen - 1, end - 1 + args.flank)
        cols = colmap[lo:hi + 1]
        sub = slice_columns(aseqs, cols)
        snames, sseqs, dropped = drop_empty_rows(anames, sub,
                                                 args.min_residues)
        if len(sseqs) < 4:
            print(f"motif {num}: only {len(sseqs)} sequences in the region, "
                  f"skipped", file=sys.stderr)
            continue

        md = os.path.join(args.out, f"motif_{num}")
        os.makedirs(md, exist_ok=True)
        row = {"Motif": num, "start": start, "end": end,
               "region_columns": len(cols), "n_seq": len(sseqs),
               "dropped_empty_rows": dropped}

        # structure, projected onto the extracted columns via the reference
        ss = None
        if args.structures:
            cand = [f for f in os.listdir(args.structures)
                    if f.endswith(".dbn")
                    and (MOTIF_RE.search(f) or (None,))
                    and MOTIF_RE.search(f)
                    and int(MOTIF_RE.search(f).group(1)) == num]
            if cand:
                raw = read_dbn(os.path.join(args.structures, cand[0]))
                if raw:
                    raw = normalise(raw)
                    ref_row = sseqs[snames.index(anames[ref_i])] \
                        if anames[ref_i] in snames else None
                    if ref_row is None:
                        print(f"motif {num}: reference absent from the "
                              f"region, cannot place the structure",
                              file=sys.stderr)
                    else:
                        rescols = [i for i, c in enumerate(ref_row)
                                   if c not in GAPS]
                        if len(rescols) != len(raw):
                            print(f"motif {num}: structure is {len(raw)} nt "
                                  f"but the reference has {len(rescols)} nt "
                                  f"in this region; skipped", file=sys.stderr)
                        else:
                            tmp = ["."] * len(ref_row)
                            for k, c in enumerate(rescols):
                                tmp[c] = raw[k]
                            ss = rebalance("".join(tmp))
                            row["base_pairs"] = len(pairs_of(ss))
            else:
                print(f"motif {num}: no .dbn found", file=sys.stderr)

        if ss:
            sseqs, ss, emptied = drop_empty_columns(sseqs, ss)
            row["dropped_empty_cols"] = emptied
            row["base_pairs_tested"] = len(pairs_of(ss))

        sto = os.path.join(md, f"motif_{num}_agnostic.sto")
        write_sto(sto, snames, sseqs, ss)
        row["alignment_file"] = os.path.relpath(sto, args.out)

        if ss:
            od = os.path.join(md, "rscape")
            run([args.rscape_bin] + args.rscape_args.split() + ["--outdir", od,
                 sto], False, os.path.join(md, "rscape.log"))
            cov, exp = read_rscape(od)
            row["cov_agnostic"] = "" if cov is None else cov
            row["expected_agnostic"] = "" if exp is None else exp
            print(f"motif {num}: {len(sseqs)} seqs, "
                  f"{row.get('base_pairs_tested')} bp tested, "
                  f"covarying = {cov}")
        else:
            print(f"motif {num}: region extracted, no structure tested")
        rows.append(row)

    if not rows:
        sys.exit("No motif could be processed.")

    # ---- 7. comparison ---------------------------------------------------
    if args.compare and os.path.isfile(args.compare):
        prev = {}
        for r in read_table(args.compare):
            low = {(k or "").strip().lower(): (v or "").strip()
                   for k, v in r.items() if k}
            num = None
            for k, v in low.items():
                if "motif" in k and v:
                    m = re.search(r"(\d+)", v)
                    if m:
                        num = int(m.group(1))
                        break
            if num is None:
                continue
            for k in ("cov_published", "cov_original", "covarying bp",
                      "covarying_bp"):
                if low.get(k):
                    try:
                        prev[num] = int(float(low[k]))
                    except ValueError:
                        pass
                    break
        hit = 0
        for r in rows:
            if r["Motif"] in prev:
                r["cov_cm_based"] = prev[r["Motif"]]
                hit += 1
                a = r.get("cov_agnostic")
                if isinstance(a, int):
                    b = prev[r["Motif"]]
                    r["change"] = ("retained" if a > 0 and b > 0 else
                                   "LOST" if b > 0 and a == 0 else
                                   "gained" if a > 0 and b == 0 else
                                   "none either way")
        print(f"\nMatched {hit} motif(s) to {args.compare}")

    cols = ["Motif", "start", "end", "region_columns", "n_seq",
            "dropped_empty_rows", "dropped_empty_cols", "base_pairs",
            "base_pairs_tested", "cov_agnostic", "expected_agnostic",
            "cov_cm_based", "change", "alignment_file"]
    tpath = os.path.join(args.out, "agnostic_covariation.tsv")
    with open(tpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {tpath}")

    cov = [r for r in rows if isinstance(r.get("cov_agnostic"), int)]
    pos = [r for r in cov if r["cov_agnostic"] > 0]
    print(f"\n{len(pos)} of {len(cov)} motifs retain significant covariation "
          f"when the alignment is built without reference to any structure.")
    if pos:
        print(f"  {', '.join(str(r['Motif']) for r in pos)}")
    lost = [r for r in rows if r.get("change") == "LOST"]
    if lost:
        print(f"\n{len(lost)} motif(s) lost covariation relative to the "
              f"CM-based analysis: "
              f"{', '.join(str(r['Motif']) for r in lost)}")
        print("  Their earlier support came from the alignment procedure.")
    depths = Counter(r["n_seq"] for r in rows)
    if len(depths) == 1:
        print(f"\nAlignment depth is now uniform at "
              f"{next(iter(depths))} sequences for every motif.")
    else:
        d = [r["n_seq"] for r in rows]
        print(f"\nAlignment depth {min(d)}-{max(d)} sequences "
              f"(varies where genomes are gapped in a region).")
    print("\nNext: recompute the ti/tv analysis on these same alignments, "
          "which are also free of the circularity:")
    print(f"  titv_table.py {args.out}/motif_*/*_agnostic.sto "
          f"--out titv_agnostic --weight")


if __name__ == "__main__":
    main()
