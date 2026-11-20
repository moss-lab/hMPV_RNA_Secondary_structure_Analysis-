#!/usr/bin/env python3
"""Earlier transition/transversion table builder. Superseded by titv_table.py."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys

MOTIF_RE = re.compile(r"motif[_-]?(\d+)(?=[._\-]|$)", re.IGNORECASE)

OUT_COLUMNS = [
    "Class", "Motif", "Gene", "Start", "End", "MFE (kcal/mol)", "z-score",
    "ED", "n_seq", "n_subs", "Stem ti/tv", "Loop ti/tv",
    "OR (stem vs loop)", "OR 95% CI", "Fisher P", "q (BH)",
    "Direction", "Flags",
]

NUMERIC_RIGHT = {
    "Motif", "Start", "End", "MFE (kcal/mol)", "z-score", "ED", "n_seq",
    "n_subs", "Stem ti/tv", "Loop ti/tv", "OR (stem vs loop)",
    "Fisher P", "q (BH)",
}

# Fallback metadata: gene, motif, coordinates, MFE, z-score, ED.
# Rounded, so --metadata with the full-precision sheet is preferred.
BUILTIN_META = """N,1,96-194,-23,-2.99,17.63
N,2,195-209,-2.1,-4.75,0.6
N,3,477-575,-18.5,-1.7,15.02
N,4,871-973,-28.8,-2.08,14.93
N,5,1018-1026,0,0.15,0.53
N,6,1045-1067,-5.6,-2.75,0.38
N,7,1082-1096,-2.6,-3.47,0.89
N,8,1109-1122,-0.1,-0.6,2.22
N,9,1134-1167,-13.1,-4.05,1.59
N,11,1198-1207,-1.3,-2.93,0.63
P,12,1212-1227,0,0.2,3.04
P,13,1251-1291,-12.4,-4.42,9.24
P,14,1577-1649,-13,-2.56,7.95
P,15,1687-1767,-15.6,-2.23,20.15
P,16,1804-1916,-20.5,-0.94,19.4
P,17,1985-2078,-10.3,-1.59,10.07
M,18,2314-2330,-0.4,-0.63,1.71
M,19,2347-2404,-16.8,-2.96,2.93
M,20,2483-2556,-17.2,-2.61,5.04
M,21,2561-2575,-4.3,-2.93,0.4
M,22,2588-2618,-8,-2.88,1.71
M,24,2660-2682,-4.5,-3.16,0.54
M,25,2691-2720,-11.8,-3.89,0.63
M,26,2738-2797,-12.1,-1.97,5.45
M,27,2809-2866,-17.8,-2.15,4.91
F,28,3009-3092,-16.2,-2.11,5.88
F,29,3103-3117,-1.3,-1.02,1.58
F,30,3370-3425,-21.6,-2.88,7.45
F,31,3596-3605,0,0.23,1.43
F,32,3610-3676,-17.5,-1.81,9.14
F,33,3682-3698,-6,-4.73,0.06
F,34,3754-3787,-9.1,-2.77,1.69
F,35,3789-3819,-4.6,-1.17,2.16
F,36,3944-3979,-9.6,-2.3,3.03
F,37,4006-4030,0,0.72,0.49
F,38,4038-4054,-3.3,-2.16,1.06
F,39,4061-4084,-10.6,-3.16,0.19
F,40,4086-4101,-4.4,-3.42,0.09
F,41,4269-4371,-16.8,-1.14,15.78
F,42,4429-4465,-12.7,-3.2,0.45
F,43,4574-4612,-2.5,-0.92,8.07
F,44,4615-4625,-2.4,-3.67,0.49
M2,45,4628-4673,-8.5,-1.18,5.51
M2,46,4724-4743,-5.1,-3.2,0.35
M2,47,4744-4781,-15.2,-3.36,0.78
M2,48,4792-4811,-1.1,-0.01,2.24
M2,49,4817-4830,-3.4,-1.99,0.21
M2,50,4918-4974,-14.1,-1.4,8.05
M2,51,4977-4995,-3.3,-2.15,0.69
M2,52,4998-5012,-0.1,0.11,1.59
M2,53,5016-5074,-9.6,-1.27,9.48
M2,54,5255-5301,-10.5,-2.07,3.73
M2,55,5426-5444,-1,-5.11,1.57
SH,56,5481-5569,-9.5,-0.56,18.34
SH,57,5576-5614,-7.9,-2.66,11.86
SH,58,5730-5745,-2.2,-4.25,0.55
SH,59,5757-5823,-6.8,-0.8,6.52
SH,60,5848-5935,-11,-2.19,9.53
SH,61,6015-6032,-1.2,-2.39,1.26
SH,62,6042-6096,-13.3,-3.0,2.97
G,63,6100-6116,0,0.24,3.33
G,64,6131-6144,-2,-7.07,0.83
G,65,6149-6169,-4.1,-4.92,0.77
G,66,6191-6254,-9.77,-2.89,6.2
L,67,6537-6639,-12.7,-0.47,17.09
L,68,7043-7068,-2.2,-0.53,1.69
L,69,7093-7139,-4.48,-3.42,3.62
L,70,7155-7200,-13.1,-4.33,3.07
L,71,7220-7239,-6.4,-4.25,1
L,72,7293-7303,0,0.13,0.24
L,73,7508-7602,-19.4,-1.75,13.97
L,74,7662-7711,-12.5,-2.62,1.85
L,75,7810-7844,-6.1,-1.78,1.79
L,76,7848-7873,-3.2,-1.33,2.43
L,77,7900-7948,-6,-1.12,10.34
L,78,7963-8003,-13.7,-3.58,1.05
L,79,8314-8390,-15.6,-1.99,8.6
L,80,8663-8749,-16.3,-2.06,8.04
L,81,9024-9047,-1.6,-1.16,2
L,82,9054-9127,-20.4,-2.67,6.46
L,83,9310-9325,0,0.4,1.98
L,84,9338-9418,-21.9,-4.06,2.92
L,85,9427-9480,-14,-2.74,6.07
L,86,9912-9925,-2.1,-3.48,0.46
L,87,9929-9982,-18.9,-4.12,7.09
L,88,9984-9998,-2.9,-2.54,0.65
L,89,10002-10027,-8.6,-3.09,1.15
L,90,10036-10061,-4.4,-2.43,1.19
L,91,10470-10553,-20.4,-2.82,2.89
L,92,10555-10563,-1.2,-5.05,1.97
L,93,10853-10886,-6.9,-2.04,3.44
L,94,10887-10947,-8,-0.53,16.01
L,95,11122-11208,-17.4,-2.85,16.29
L,96,11261-11344,-19,-2.48,6.16
L,97,11351-11361,-0.7,-2.77,1.39
L,98,11785-11863,-19.9,-1.99,11.65
L,99,11867-11878,0,0.15,1.08
L,100,12071-12123,-6,-0.87,3.75
L,101,12128-12213,-22,-3.08,7.69
L,102,12236-12246,0,0.18,0.95
L,103,12294-12392,-20.8,-2.23,14.38
L,104,12842-12871,-6.7,-2.4,2.33
L,105,12905-12946,-7.9,-2.66,3.35
L,106,13197-13309,-14.3,-1.66,21.04"""


def to_float(v):
    try:
        f = float(str(v).strip())
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def motif_number(text: str):
    """Motif number from 'nc_hmnv_crna_motif_3', 'Motif 3', or a bare '3'."""
    text = str(text).strip()
    if not text:
        return None
    m = MOTIF_RE.search(text)
    if m:
        return int(m.group(1))
    if re.fullmatch(r"\d+", text):
        return int(text)
    # a lone number inside something like "Motif no. 3" or "3 (N)"
    m = re.fullmatch(r"[^0-9]*?(\d+)[^0-9]*", text)
    return int(m.group(1)) if m else None


def loose(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def sniff_read(path: str):
    """Rows of a CSV or TSV as dicts with lowercased, stripped keys."""
    with open(path, newline="", errors="replace") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        for row in csv.DictReader(fh, delimiter=delim):
            yield {(k or "").strip().lower(): (v or "").strip()
                   for k, v in row.items() if k is not None}


def pick(row: dict, *names, contains=False):
    for n in names:
        if n in row and row[n] != "":
            return row[n]
    if contains:
        for key, val in row.items():
            if val != "" and any(n in key for n in names):
                return val
    return ""


def load_metadata(path):
    """{motif number: {gene, start, end, mfe, z, ed, cls}}"""
    meta = {}
    if path:
        for row in sniff_read(path):
            num = motif_number(pick(row, "motif", "motif_n", "motif no",
                                    "motif no.", "structure", contains=True))
            if num is None:
                continue
            start = pick(row, "start")
            end = pick(row, "end")
            if not start or not end:
                coords = pick(row, "coordinates", "coordinate", contains=True)
                m = re.search(r"(\d+)\s*[-–]\s*(\d+)", coords)
                if m:
                    start, end = m.group(1), m.group(2)
            meta[num] = {
                "gene": pick(row, "gene"),
                "start": start, "end": end,
                "mfe": pick(row, "mfe", "mfe_kcal_mol", "mfe (kcal/mol)",
                            contains=True),
                "z": pick(row, "z_score", "z-score", "zscore", "z_score_table",
                          contains=True),
                "ed": pick(row, "ed"),
                "cls": pick(row, "class", "classo", "tier", contains=True),
            }
        print(f"Loaded metadata for {len(meta)} motifs from {path}")
    else:
        for line in BUILTIN_META.splitlines():
            gene, num, coords, mfe, z, ed = line.split(",")
            start, end = coords.split("-")
            meta[int(num)] = {"gene": gene, "start": start, "end": end,
                              "mfe": mfe, "z": z, "ed": ed, "cls": ""}
        print(f"Using the built-in metadata for {len(meta)} motifs. Its MFE "
              f"and z-score are rounded; pass --metadata with your full "
              f"ScanFold sheet for publication values.")
    return meta


def fmt(v, places, na="NA"):
    f = to_float(v)
    return na if f is None else f"{f:.{places}f}"


def fmt_p(v, na="NA"):
    f = to_float(v)
    if f is None:
        return na
    if f == 0:
        return "<1e-300"
    if f < 0.001:
        return f"{f:.1e}"
    return f"{f:.3f}"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("summary", help="titv_summary.csv from titv_analysis.py")
    ap.add_argument("--metadata",
                    help="CSV/TSV of motif gene/coordinates/MFE/z-score/ED")
    ap.add_argument("--out", default="table4_titv.csv")
    ap.add_argument("--markdown", help="also write an aligned markdown table")
    ap.add_argument("--q-cut", type=float, default=0.05)
    ap.add_argument("--min-subs", type=int, default=20,
                    help="below this a row is flagged low-n (default 20)")
    ap.add_argument("--clean-only", action="store_true",
                    help="keep only rows with no flags")
    ap.add_argument("--sort", choices=["motif", "q", "or"], default="motif")
    ap.add_argument("--na", default="NA")
    args = ap.parse_args()

    if not os.path.isfile(args.summary):
        sys.exit(f"No such file: {args.summary}")
    meta = load_metadata(args.metadata)

    rows, unmatched = [], []
    for rec in sniff_read(args.summary):
        structure = pick(rec, "structure", "name", contains=True)
        num = motif_number(structure)
        if num is None:
            unmatched.append(structure)
            continue
        m = meta.get(num, {})

        orv = to_float(pick(rec, "oddsratio_stem_vs_loop", "oddsratio",
                            contains=True))
        lo = to_float(pick(rec, "or_ci_low"))
        hi = to_float(pick(rec, "or_ci_high"))
        p = to_float(pick(rec, "fisher_p"))
        q = to_float(pick(rec, "q"))
        nsub = to_float(pick(rec, "n_substitutions"))
        zero_cell = str(pick(rec, "or_zero_cell")).lower() in ("true", "1",
                                                              "yes")

        flags = []
        if zero_cell:
            flags.append("zero-cell")
        if lo is not None and hi is not None and lo <= 1 <= hi:
            flags.append("CI-spans-1")
        if nsub is not None and nsub < args.min_subs:
            flags.append("low-n")
        if q is not None and q >= args.q_cut:
            flags.append("ns")

        direction = pick(rec, "direction")
        if not direction and orv is not None:
            direction = "stem-biased" if orv > 1 else "loop-biased"

        row = {
            "Class": m.get("cls", "") or pick(rec, "tier", "tier_label",
                                              contains=True),
            "Motif": num,
            "Gene": m.get("gene", ""),
            "Start": m.get("start", ""),
            "End": m.get("end", ""),
            "MFE (kcal/mol)": m.get("mfe", ""),
            "z-score": m.get("z", ""),
            "ED": m.get("ed", ""),
            "n_seq": pick(rec, "n_seq"),
            "n_subs": "" if nsub is None else int(round(nsub)),
            "Stem ti/tv": fmt(pick(rec, "stem_ti_tv"), 2, args.na),
            "Loop ti/tv": fmt(pick(rec, "loop_ti_tv"), 2, args.na),
            "OR (stem vs loop)": fmt(orv, 2, args.na),
            "OR 95% CI": (f"{lo:.2f}-{hi:.2f}"
                          if lo is not None and hi is not None else args.na),
            "Fisher P": fmt_p(p, args.na),
            "q (BH)": fmt_p(q, args.na),
            "Direction": direction,
            "Flags": "; ".join(flags) if flags else "",
        }
        for c in ("MFE (kcal/mol)", "z-score", "ED"):
            if row[c] == "":
                row[c] = args.na
        rows.append(row)

    if unmatched:
        print(f"  ! {len(unmatched)} summary rows had no motif number and "
              f"were skipped: {', '.join(unmatched[:5])}", file=sys.stderr)
    if not rows:
        sys.exit("No rows to write.")

    missing_meta = [r["Motif"] for r in rows if not r["Gene"]]
    if missing_meta:
        print(f"  ! no metadata for motifs: "
              f"{', '.join(map(str, missing_meta))}", file=sys.stderr)

    flagged = [r for r in rows if r["Flags"]]
    if args.clean_only:
        rows = [r for r in rows if not r["Flags"]]

    keys = {"motif": lambda r: r["Motif"],
            "q": lambda r: (to_float(r["q (BH)"].lstrip("<")) or 9e99),
            "or": lambda r: -(to_float(r["OR (stem vs loop)"]) or -9e99)}
    rows.sort(key=keys[args.sort])

    delim = "\t" if args.out.lower().endswith((".tsv", ".txt")) else ","
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter=delim,
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {args.out}")

    if args.markdown:
        widths = [max(len(c), *(len(str(r[c])) for r in rows))
                  for c in OUT_COLUMNS]
        def cell(c, text, w):
            return str(text).rjust(w) if c in NUMERIC_RIGHT \
                else str(text).ljust(w)
        with open(args.markdown, "w") as fh:
            fh.write("| " + " | ".join(cell(c, c, w) for c, w
                                       in zip(OUT_COLUMNS, widths)) + " |\n")
            fh.write("|" + "|".join(("-" * (w + 1) + ":") if c in NUMERIC_RIGHT
                                    else ("-" * (w + 2))
                                    for c, w in zip(OUT_COLUMNS, widths))
                     + "|\n")
            for r in rows:
                fh.write("| " + " | ".join(cell(c, r[c], w) for c, w
                                           in zip(OUT_COLUMNS, widths))
                         + " |\n")
        print(f"Wrote markdown to {args.markdown}")

    sig = [r for r in rows if not r["Flags"]]
    print(f"\n{len(sig)} of {len(rows)} rows carry no flags "
          f"(q < {args.q_cut}, CI excludes 1, >= {args.min_subs} "
          f"substitutions, no empty cell).")
    if flagged and not args.clean_only:
        counts = {}
        for r in flagged:
            for f in r["Flags"].split("; "):
                counts[f] = counts.get(f, 0) + 1
        print("Flagged rows:")
        for f, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            print(f"    {f:12s} {n}")
        spans = [r["Motif"] for r in flagged if "CI-spans-1" in r["Flags"]]
        zeros = [r["Motif"] for r in flagged if "zero-cell" in r["Flags"]]
        if zeros:
            print(f"\n  Motifs {', '.join(map(str, zeros))} have an empty "
                  f"cell in the 2x2 table. Their odds ratio comes from the "
                  f"+0.5 correction and should not be quoted as an effect "
                  f"size.")
        if spans:
            print(f"  Motifs {', '.join(map(str, spans))} have a CI that "
                  f"includes 1, so the direction is not established for "
                  f"them individually.")
    loop = [r["Motif"] for r in rows if r["Direction"] == "loop-biased"]
    if loop:
        print(f"\n  {len(loop)} structure(s) are loop-biased "
              f"({', '.join(map(str, loop))}). Keep them in the table: "
              f"dropping rows that run against the hypothesis is what "
              f"makes a table look selective.")


if __name__ == "__main__":
    main()
