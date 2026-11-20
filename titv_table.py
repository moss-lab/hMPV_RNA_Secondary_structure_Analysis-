#!/usr/bin/env python3
"""Stem versus loop transition/transversion bias with Fisher exact test and multiple-testing correction."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from collections import Counter

# Optional accelerators. The script is correct without them.
try:
    import numpy as _np
except ImportError:
    _np = None
try:
    from scipy.stats import fisher_exact as _scipy_fisher
except ImportError:
    _scipy_fisher = None

RNA = "ACGU"
TRANSITION = {"A": "G", "G": "A", "C": "U", "U": "C"}
GAPS = set("-.~_")
PAIRED_CH = set("()[]{}<>AaBbCcDd")      # WUSS / VARNA pairing symbols
LOOP_CH = set(".,:_-~")
ALN_EXT = (".sto", ".stk", ".stockholm", ".fa", ".fasta", ".fna", ".afa",
           ".aln", ".afasta")
FOURFOLD = {"GC", "CG", "GG", "CU", "CC", "UC", "AC", "GU"}


# Reading

def sniff_format(path):
    """FASTA or Stockholm, decided by content. Extensions lie."""
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
            names.append(line[1:].split()[0])
            cur = []
        elif line:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return names, [s.upper().replace("T", "U") for s in seqs], {}


def read_stockholm(path):
    """Handles interleaved blocks: a repeated name continues one sequence."""
    names, seqs, gc = [], [], {}
    index = {}
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
        p = s.rsplit(None, 1)          # names may contain spaces
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
    return (read_fasta(path) if sniff_format(path) == "fasta"
            else read_stockholm(path))


def read_dbn(path):
    """Last line that looks like dot-bracket, ignoring headers and sequence."""
    text = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    for line in reversed(text):
        s = line.strip().split()[0] if line.strip() else ""
        if not s or s.startswith(">"):
            continue
        if set(s) <= (PAIRED_CH | LOOP_CH) and any(c in "()[]{}<>" for c in s):
            return s
    raise ValueError(f"no dot-bracket line found in {path}")


# Pairing alignments with structures

def strip_known_suffixes(basename, extra_patterns=()):
    b = basename
    for _ in range(6):
        b = re.sub(r"\.(sto|stk|stockholm|fa|fasta|fna|afa|afasta|aln|dbn|txt)$",
                   "", b, flags=re.I)
        b = re.sub(r"\.(named|cacofold|r2r|sorted|stockholm(_\d+)?)$", "", b,
                   flags=re.I)
    b = b.lower().replace("-", "_")
    for pat in extra_patterns:
        b = re.sub(pat, "", b, flags=re.I)
    return re.sub(r"[._]+$", "", b)


def pair_key(path, extra_patterns=()):
    """A key two files must share to describe the same structure. Generic: it."""
    return strip_known_suffixes(os.path.basename(path), extra_patterns)


def map_structure_to_columns(dbn, seqs, gc, names, ref_name=None):
    """Project an ungapped dot-bracket onto gapped alignment columns."""
    L = len(seqs[0])
    if len(dbn) == L:
        return dbn, [True] * L, "exact (structure spans every column)"

    def place(row, label):
        cols = [i for i, c in enumerate(row) if c not in GAPS]
        if len(cols) != len(dbn):
            return None
        out = ["."] * L
        mask = [False] * L
        for k, i in enumerate(cols):
            out[i] = dbn[k]
            mask[i] = True
        return "".join(out), mask, label

    if ref_name:
        for nm, sq in zip(names, seqs):
            if ref_name.lower() in nm.lower():
                r = place(sq, f"reference sequence {nm!r}")
                if r:
                    return r
    rf = gc.get("RF")
    if rf and len(rf) == L:
        r = place(rf, "#=GC RF consensus")
        if r:
            return r
    for idx, sq in enumerate(seqs):
        r = place(sq, f"sequence {names[idx]!r}")
        if r:
            return r
    return None, None, None


# Statistics, implemented so scipy is optional

def _lbinom(n, k):
    if k < 0 or k > n:
        return -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def fisher_exact_two_sided(a, b, c, d):
    """Two-sided Fisher exact P for [[a,b],[c,d]], by summing the."""
    a, b, c, d = int(round(a)), int(round(b)), int(round(c)), int(round(d))
    if min(a, b, c, d) < 0:
        return float("nan")
    n = a + b + c + d
    if n == 0:
        return float("nan")
    r1, r2 = a + b, c + d
    c1, c2 = a + c, b + d
    if min(r1, r2, c1, c2) == 0:
        return 1.0

    log_denom = _lbinom(n, c1)

    def logp(x):
        return _lbinom(r1, x) + _lbinom(r2, c1 - x) - log_denom

    lo, hi = max(0, c1 - r2), min(c1, r1)
    obs = logp(a)
    tol = 1e-7
    total = 0.0
    for x in range(lo, hi + 1):
        lp = logp(x)
        if lp <= obs + tol:
            total += math.exp(lp)
    return min(1.0, total)


def odds_ratio(a, b, c, d, correction=0.5):
    """Haldane-Anscombe corrected odds ratio, so an empty cell gives a finite."""
    A, B, C, D = a + correction, b + correction, c + correction, d + correction
    if B * C == 0:
        return float("nan")
    return (A * D) / (B * C)


def or_ci_woolf(a, b, c, d, correction=0.5, z=1.959963984540054):
    """Woolf log-OR interval. Analytic, always available."""
    A, B, C, D = a + correction, b + correction, c + correction, d + correction
    if min(A, B, C, D) <= 0:
        return float("nan"), float("nan")
    lor = math.log((A * D) / (B * C))
    se = math.sqrt(1 / A + 1 / B + 1 / C + 1 / D)
    return math.exp(lor - z * se), math.exp(lor + z * se)


def or_ci_bootstrap(table, n_boot=5000, seed=0, correction=0.5):
    """Percentile bootstrap on the odds ratio, resampling the 2x2 multinomial."""
    a, b, c, d = table
    n = int(round(a + b + c + d))
    if _np is None or n <= 0:
        return or_ci_woolf(a, b, c, d, correction) + ("Woolf",)
    rng = _np.random.default_rng(seed)
    flat = _np.array([a, b, c, d], float)
    if flat.sum() <= 0:
        return float("nan"), float("nan"), "none"
    draws = rng.multinomial(n, flat / flat.sum(), size=n_boot).astype(float)
    draws += correction
    ors = (draws[:, 0] * draws[:, 3]) / (draws[:, 1] * draws[:, 2])
    ors = ors[_np.isfinite(ors)]
    if ors.size < 100:
        return or_ci_woolf(a, b, c, d, correction) + ("Woolf",)
    lo, hi = _np.percentile(ors, [2.5, 97.5])
    return float(lo), float(hi), f"bootstrap n={n_boot}"


def benjamini_hochberg(pvals):
    """FDR q values, NaN-safe, order preserved."""
    out = [float("nan")] * len(pvals)
    idx = [i for i, p in enumerate(pvals)
           if p is not None and isinstance(p, float) and math.isfinite(p)]
    if not idx:
        return out
    m = len(idx)
    ordered = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank, i in enumerate(reversed(ordered)):
        q = pvals[i] * m / (m - rank)
        prev = min(prev, q)
        out[i] = min(1.0, prev)
    return out


# Counting

def henikoff_weights(seqs):
    """Position-based weights; damps phylogenetic redundancy."""
    n, L = len(seqs), len(seqs[0])
    w = [0.0] * n
    for i in range(L):
        col = [s[i] for s in seqs]
        obs = [c for c in col if c in RNA]
        if not obs:
            continue
        tally = Counter(obs)
        r = len(tally)
        if r < 2:
            continue
        for k, c in enumerate(col):
            if c in RNA:
                w[k] += 1.0 / (r * tally[c])
    total = sum(w)
    if total <= 0:
        return [1.0] * n
    return [x * n / total for x in w]


def consensus_of(seqs, weights):
    L = len(seqs[0])
    out = []
    for i in range(L):
        tally = Counter()
        for s, ww in zip(seqs, weights):
            if s[i] in RNA:
                tally[s[i]] += ww
        out.append(tally.most_common(1)[0][0] if tally else "-")
    return "".join(out)


def classify_columns(dbn):
    """1 paired, 0 unpaired, -1 ignore."""
    return [1 if ch in PAIRED_CH else (0 if ch in LOOP_CH else -1)
            for ch in dbn]


def codon_positions(L, frame):
    return [((i + frame) % 3) + 1 for i in range(L)]


def fourfold_mask(consensus, frame):
    L = len(consensus)
    cp = codon_positions(L, frame)
    mask = [False] * L
    for i in range(L):
        if cp[i] == 3 and i >= 2 and consensus[i - 2:i] in FOURFOLD:
            mask[i] = True
    return mask


def count_substitutions(seqs, consensus, cls, weights, use_cols):
    """Per-mutation-type counts plus the pooled 2x2 table."""
    rows = []
    pooled = {"stem": {"ti": 0.0, "tv": 0.0, "opp": 0.0},
              "loop": {"ti": 0.0, "tv": 0.0, "opp": 0.0}}
    L = len(consensus)

    # Index columns once per (class, consensus base) to keep this O(L*N).
    buckets = {(w, nt): [] for w in (0, 1) for nt in RNA}
    for i in range(L):
        if use_cols[i] and cls[i] in (0, 1) and consensus[i] in RNA:
            buckets[(cls[i], consensus[i])].append(i)

    for nt in RNA:
        for tgt in RNA:
            if nt == tgt:
                continue
            kind = "Transition" if TRANSITION[nt] == tgt else "Transversion"
            rec = {"Mutation": f"{nt}->{tgt}", "From": nt, "Type": kind}
            for label, want in (("stem", 1), ("loop", 0)):
                m = opp = 0.0
                for i in buckets[(want, nt)]:
                    for s, ww in zip(seqs, weights):
                        ch = s[i]
                        if ch not in RNA:
                            continue
                        opp += ww
                        if ch == tgt:
                            m += ww
                rec[f"count_{label}"] = m
                rec[f"opportunity_{label}"] = opp
                rec[f"rate_{label}"] = m / opp if opp else float("nan")
            rows.append(rec)

    for key in ("stem", "loop"):
        ti = sum(r[f"count_{key}"] for r in rows if r["Type"] == "Transition")
        tv = sum(r[f"count_{key}"] for r in rows
                 if r["Type"] == "Transversion")
        # Opportunity is recorded once per (from, to) pair; the four
        # transition rows cover each source base exactly once, so their sum
        # is the per-site opportunity.
        opp = sum(r[f"opportunity_{key}"] for r in rows
                  if r["Type"] == "Transition")
        pooled[key] = {"ti": ti, "tv": tv, "opp": opp}
    return rows, pooled


def conservation(seqs, consensus, cls, weights, use_cols):
    st, lp = [], []
    W = sum(weights)
    for i in range(len(consensus)):
        if not use_cols[i] or cls[i] < 0:
            continue
        match = sum(w for s, w in zip(seqs, weights) if s[i] == consensus[i])
        (st if cls[i] == 1 else lp).append(match / W if W else float("nan"))
    sc = sum(st) / len(st) if st else float("nan")
    lc = sum(lp) / len(lp) if lp else float("nan")
    return {"stem_conservation": sc, "loop_conservation": lc,
            "n_stem_cols": len(st), "n_loop_cols": len(lp)}


# One structure

def analyse(aln_path, dbn_path, args):
    name = pair_key(aln_path, args._strip)
    names, seqs, gc = read_alignment(aln_path)
    if not seqs:
        return None, f"{name}: no sequences"

    L = len(seqs[0])
    ragged = [(n, len(s)) for n, s in zip(names, seqs) if len(s) != L]
    if ragged:
        lens = sorted({len(s) for s in seqs})
        return None, (f"{name}: sequences are not all the same length "
                      f"({len(seqs)} seqs, lengths {lens[:6]}), parsed as "
                      f"{sniff_format(aln_path).upper()}")

    dbn = read_dbn(dbn_path) if dbn_path else gc.get("SS_cons")
    if not dbn:
        return None, f"{name}: no .dbn supplied and no #=GC SS_cons present"

    mapped, ref_mask, how = map_structure_to_columns(
        dbn, seqs, gc, names, args.ref)
    if mapped is None:
        ung = sorted({sum(1 for c in s if c not in GAPS) for s in seqs})
        return None, (f"{name}: structure is {len(dbn)} nt but the alignment "
                      f"has {L} columns with ungapped lengths {ung[:6]}; no "
                      f"row matches, cannot map")
    dbn = mapped

    weights = henikoff_weights(seqs) if args.weight else [1.0] * len(seqs)
    consensus = consensus_of(seqs, weights)
    cls = classify_columns(dbn)

    frame = args._frames.get(re.sub(r"[^a-z0-9]", "", name.lower()),
                             args.frame)
    use = list(ref_mask)
    if args.restrict == "third":
        cp = codon_positions(L, frame)
        use = [u and cp[i] == 3 for i, u in enumerate(use)]
    elif args.restrict == "fourfold":
        ff = fourfold_mask(consensus, frame)
        use = [u and ff[i] for i, u in enumerate(use)]
    use = [u and consensus[i] in RNA for i, u in enumerate(use)]
    if not any(use):
        return None, f"{name}: no usable columns after --restrict {args.restrict}"

    mut, pooled = count_substitutions(seqs, consensus, cls, weights, use)
    s, l = pooled["stem"], pooled["loop"]
    a, b, c, d = s["ti"], s["tv"], l["ti"], l["tv"]

    zero_cell = min(a, b, c, d) == 0
    p = fisher_exact_two_sided(a, b, c, d)
    orv = odds_ratio(a, b, c, d)
    lo, hi, ci_method = or_ci_bootstrap((a, b, c, d), args.boot, args.seed)

    row = {
        "Structure": name,
        "alignment": os.path.basename(aln_path),
        "structure_file": os.path.basename(dbn_path) if dbn_path
                          else "SS_cons",
        "mapping": how,
        "n_seq": len(seqs),
        "n_cols": L,
        "cols_used": sum(1 for u in use if u),
        "restrict": args.restrict,
        "weighted": bool(args.weight),
        "stem_ti": a, "stem_tv": b, "loop_ti": c, "loop_tv": d,
        "n_subs": a + b + c + d,
        "stem_opportunity": s["opp"], "loop_opportunity": l["opp"],
        "stem_ti_tv": a / b if b else float("nan"),
        "loop_ti_tv": c / d if d else float("nan"),
        # neutral expectation: 1 transition and 2 transversions per site,
        # so a normalised value of 1.0 means no bias
        "stem_ti_tv_norm": 2 * a / b if b else float("nan"),
        "loop_ti_tv_norm": 2 * c / d if d else float("nan"),
        "stem_sub_rate": (a + b) / (3 * s["opp"]) if s["opp"] else float("nan"),
        "loop_sub_rate": (c + d) / (3 * l["opp"]) if l["opp"] else float("nan"),
        "odds_ratio": orv,
        "or_ci_low": lo, "or_ci_high": hi, "ci_method": ci_method,
        "fisher_p": p,
        "zero_cell": zero_cell,
        "direction": ("stem-biased" if orv == orv and orv > 1
                      else "loop-biased" if orv == orv else "undefined"),
        **conservation(seqs, consensus, cls, weights, use),
    }
    return (row, mut), None


# Metadata and output

def read_table(path):
    with open(path, newline="", errors="replace") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        return list(csv.DictReader(fh, delimiter=delim))


def join_key(text, mode):
    """How a metadata row is matched to a structure. 'auto' pulls the last."""
    t = str(text).strip()
    if mode == "name":
        return re.sub(r"[^a-z0-9]", "", t.lower())
    nums = re.findall(r"(\d+)", t)
    return nums[-1] if nums else re.sub(r"[^a-z0-9]", "", t.lower())


def fmt_num(v, places):
    if v is None:
        return ""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(f):
        return "NA"
    if places is None:
        return str(int(f)) if f == int(f) else f"{f:.1f}"
    return f"{f:.{places}f}"


def fmt_p(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "NA"
    if not math.isfinite(f):
        return "NA"
    if f == 0:
        return "<1e-300"
    return f"{f:.1e}" if f < 0.001 else f"{f:.3f}"


DECIMALS = {"stem_ti_tv": 2, "loop_ti_tv": 2, "stem_ti_tv_norm": 2,
            "loop_ti_tv_norm": 2, "odds_ratio": 2, "or_ci_low": 2,
            "or_ci_high": 2, "stem_conservation": 3, "loop_conservation": 3,
            "stem_sub_rate": 4, "loop_sub_rate": 4,
            # counts: auto, so whole numbers print as integers and
            # weighted (fractional) counts keep their decimals
            "stem_ti": None, "stem_tv": None, "loop_ti": None,
            "loop_tv": None, "n_subs": None}

DEFAULT_COLUMNS = ["Structure", "n_seq", "n_subs", "stem_ti_tv",
                   "loop_ti_tv", "odds_ratio", "or_95_ci", "fisher_p", "q",
                   "direction", "flags"]


def write_table(rows, columns, path, markdown=None, na="NA"):
    def cell(row, col):
        if col == "or_95_ci":
            lo, hi = row.get("or_ci_low"), row.get("or_ci_high")
            try:
                return (f"{float(lo):.2f}-{float(hi):.2f}"
                        if math.isfinite(float(lo)) and math.isfinite(float(hi))
                        else na)
            except (TypeError, ValueError):
                return na
        v = row.get(col, "")
        if col in ("fisher_p", "q"):
            return fmt_p(v)
        if isinstance(v, bool):
            return "yes" if v else "no"
        if isinstance(v, (int, float)):
            return fmt_num(v, DECIMALS.get(col))
        if col == "flags":
            return str(v)          # blank means nothing is wrong with the row
        return na if v in (None, "") else str(v)

    delim = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=delim, lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow([cell(r, c) for c in columns])
    print(f"Wrote {len(rows)} rows to {path}")

    if markdown:
        text = [[cell(r, c) for c in columns] for r in rows]
        widths = [max(len(columns[i]), *(len(t[i]) for t in text))
                  if text else len(columns[i]) for i in range(len(columns))]
        numeric = {i for i in range(len(columns))
                   if all(re.fullmatch(r"[-+0-9.eE<>]+", t[i] or "x")
                          for t in text)} if text else set()
        def pad(i, s):
            return s.rjust(widths[i]) if i in numeric else s.ljust(widths[i])
        with open(markdown, "w") as fh:
            fh.write("| " + " | ".join(pad(i, c) for i, c
                                       in enumerate(columns)) + " |\n")
            fh.write("|" + "|".join(("-" * (widths[i] + 1) + ":")
                                    if i in numeric else
                                    ("-" * (widths[i] + 2))
                                    for i in range(len(columns))) + "|\n")
            for t in text:
                fh.write("| " + " | ".join(pad(i, t[i]) for i
                                           in range(len(columns))) + " |\n")
        print(f"Wrote markdown to {markdown}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="alignment files (.sto/.fa/...); .dbn files may be "
                         "mixed in and will be paired by name")
    ap.add_argument("--dbn-dir", help="directory of .dbn files")
    ap.add_argument("--dbn", help="one .dbn applied to a single alignment")
    ap.add_argument("--out", default="titv_results",
                    help="output directory (default titv_results)")
    ap.add_argument("--table", default="titv_table.csv",
                    help="final table filename inside --out")
    ap.add_argument("--markdown", action="store_true",
                    help="also write the table as aligned markdown")
    ap.add_argument("--per-structure", action="store_true",
                    help="also write a per-mutation-type CSV per structure")

    ap.add_argument("--restrict", choices=["all", "third", "fourfold"],
                    default="all",
                    help="all columns (default), third codon positions, or "
                         "fourfold-degenerate sites. Use one of the latter "
                         "for coding regions to remove protein-level "
                         "constraint.")
    ap.add_argument("--frame", type=int, default=0, choices=[0, 1, 2])
    ap.add_argument("--frames", help="CSV with columns name,frame")
    ap.add_argument("--ref", help="substring of the sequence the structure "
                                  "describes, for mapping onto columns")
    ap.add_argument("--weight", action="store_true",
                    help="Henikoff sequence weighting")
    ap.add_argument("--boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--strip", action="append", default=[],
                    help="extra regex stripped from filenames when pairing "
                         "an alignment with its structure. Repeatable.")

    ap.add_argument("--metadata", help="CSV/TSV joined into the final table")
    ap.add_argument("--meta-key", default="Motif",
                    help="column in --metadata holding the join key")
    ap.add_argument("--meta-match", choices=["auto", "name"], default="auto",
                    help="auto matches on the trailing number in each name "
                         "(motif_12 <-> 12); name matches the whole string")
    ap.add_argument("--meta-columns",
                    help="comma-separated metadata columns to place at the "
                         "front of the table (default: all of them)")
    ap.add_argument("--columns",
                    help="comma-separated list overriding the output columns "
                         "entirely. Any field computed by the script is "
                         "available; --list-fields prints them.")
    ap.add_argument("--list-fields", action="store_true",
                    help="print every available column name and exit")

    ap.add_argument("--q-cut", type=float, default=0.05)
    ap.add_argument("--min-subs", type=int, default=20,
                    help="rows with fewer substitutions are flagged low-n")
    ap.add_argument("--sort", choices=["input", "key", "q", "or"],
                    default="key")
    args = ap.parse_args()

    if args.list_fields:
        print("Per-structure fields:\n  " + "\n  ".join([
            "Structure", "alignment", "structure_file", "mapping", "n_seq",
            "n_cols", "cols_used", "restrict", "weighted", "stem_ti",
            "stem_tv", "loop_ti", "loop_tv", "n_subs", "stem_opportunity",
            "loop_opportunity", "stem_ti_tv", "loop_ti_tv",
            "stem_ti_tv_norm", "loop_ti_tv_norm", "stem_sub_rate",
            "loop_sub_rate", "odds_ratio", "or_ci_low", "or_ci_high",
            "or_95_ci", "ci_method", "fisher_p", "q", "zero_cell",
            "direction", "flags", "stem_conservation", "loop_conservation",
            "n_stem_cols", "n_loop_cols"]))
        return

    args._strip = args.strip
    args._frames = {}
    if args.frames and os.path.isfile(args.frames):
        for row in read_table(args.frames):
            low = {k.lower(): v for k, v in row.items() if k}
            nm = low.get("name") or low.get("structure") or low.get("id")
            fr = low.get("frame") or low.get("offset")
            if nm and fr not in (None, ""):
                args._frames[re.sub(r"[^a-z0-9]", "", nm.lower())] = int(fr) % 3
        print(f"Loaded reading frames for {len(args._frames)} structures")

    # ---- pair alignments with structures -------------------------------
    pool = list(args.inputs)
    if args.dbn:
        pool.append(args.dbn)
    if args.dbn_dir:
        pool += glob.glob(os.path.join(args.dbn_dir, "*.dbn"))
    alns = [f for f in pool if f.lower().endswith(ALN_EXT)]
    dbns = [f for f in pool if f.lower().endswith(".dbn")]
    if not alns:
        sys.exit("No alignment files given. Pass .sto/.fa files; .dbn alone "
                 "cannot be analysed.")

    dmap, dupes = {}, []
    for d in dbns:
        k = pair_key(d, args._strip)
        if k in dmap:
            dupes.append((k, dmap[k], d))
        else:
            dmap[k] = d
    jobs = [(a, dmap.get(pair_key(a, args._strip))) for a in alns]

    print(f"{len(alns)} alignment(s), {len(dmap)} structure file(s), "
          f"restrict={args.restrict}, "
          f"weighting={'on' if args.weight else 'off'}")
    if dupes:
        print(f"  ! {len(dupes)} duplicate structure keys; the first was kept:",
              file=sys.stderr)
        for k, f1, f2 in dupes[:5]:
            print(f"      {k}: kept {os.path.basename(f1)}, ignored "
                  f"{os.path.basename(f2)}", file=sys.stderr)
    nomatch = [a for a, d in jobs if d is None]
    if nomatch and dbns:
        print(f"  ! {len(nomatch)} alignment(s) had no matching .dbn and will "
              f"fall back to SS_cons:", file=sys.stderr)
        for a in nomatch[:5]:
            print(f"      {os.path.basename(a)} needs key "
                  f"{pair_key(a, args._strip)!r}", file=sys.stderr)

    # ---- analyse --------------------------------------------------------
    os.makedirs(args.out, exist_ok=True)
    rows, muts, problems = [], {}, []
    for aln, dbn in jobs:
        try:
            result, err = analyse(aln, dbn, args)
        except Exception as exc:
            result, err = None, f"{os.path.basename(aln)}: {exc}"
        if err:
            problems.append(err)
            print(f"  ! {err}", file=sys.stderr)
            continue
        row, mut = result
        rows.append(row)
        muts[row["Structure"]] = mut
        star = ("***" if row["fisher_p"] < 1e-3 else
                "**" if row["fisher_p"] < 1e-2 else
                "*" if row["fisher_p"] < 0.05 else "ns") \
            if math.isfinite(row["fisher_p"]) else "NA"
        print(f"  {row['Structure'][:38]:<40} "
              f"stem {fmt_num(row['stem_ti_tv'], 2):>6} "
              f"loop {fmt_num(row['loop_ti_tv'], 2):>6} "
              f"OR {fmt_num(row['odds_ratio'], 2):>7}  "
              f"P {fmt_p(row['fisher_p']):>9} {star}")

    if not rows:
        sys.exit("\nNothing could be analysed.")

    # ---- FDR and flags --------------------------------------------------
    for r, q in zip(rows, benjamini_hochberg([r["fisher_p"] for r in rows])):
        r["q"] = q
        flags = []
        if r["zero_cell"]:
            flags.append("zero-cell")
        lo, hi = r["or_ci_low"], r["or_ci_high"]
        if (isinstance(lo, float) and isinstance(hi, float)
                and math.isfinite(lo) and math.isfinite(hi)
                and lo <= 1 <= hi):
            flags.append("CI-spans-1")
        if r["n_subs"] < args.min_subs:
            flags.append("low-n")
        if not (isinstance(q, float) and math.isfinite(q) and q < args.q_cut):
            flags.append("ns")
        r["flags"] = "; ".join(flags)

    # ---- metadata join --------------------------------------------------
    meta_cols = []
    if args.metadata:
        if not os.path.isfile(args.metadata):
            sys.exit(f"No such file: {args.metadata}")
        meta_rows = read_table(args.metadata)
        if not meta_rows:
            sys.exit(f"{args.metadata} has no rows")
        headers = [h for h in meta_rows[0] if h]
        key_col = next((h for h in headers
                        if h.strip().lower() == args.meta_key.strip().lower()),
                       None)
        if key_col is None:
            sys.exit(f"--meta-key {args.meta_key!r} is not a column in "
                     f"{args.metadata}. Columns: {', '.join(headers)}")
        wanted = ([c.strip() for c in args.meta_columns.split(",")]
                  if args.meta_columns else
                  [h for h in headers if h != key_col])
        unknown = [c for c in wanted if c not in headers]
        if unknown:
            sys.exit(f"--meta-columns not in {args.metadata}: "
                     f"{', '.join(unknown)}. Columns: {', '.join(headers)}")
        index = {join_key(m[key_col], args.meta_match): m for m in meta_rows}
        hit = 0
        for r in rows:
            m = index.get(join_key(r["Structure"], args.meta_match))
            if m:
                hit += 1
                for c in wanted:
                    r[c] = m.get(c, "")
            else:
                for c in wanted:
                    r.setdefault(c, "")
        meta_cols = wanted
        print(f"\nMetadata: matched {hit} of {len(rows)} structures on "
              f"{key_col!r} ({args.meta_match} matching)")
        if hit == 0:
            print("  ! nothing matched. Check --meta-key and --meta-match; "
                  "with 'auto' both sides match on their trailing number.",
                  file=sys.stderr)

    # ---- order and write ------------------------------------------------
    if args.sort == "key":
        def sort_key(r):
            nums = re.findall(r"(\d+)", r["Structure"])
            return (int(nums[-1]) if nums else 0, r["Structure"])
        rows.sort(key=sort_key)
    elif args.sort == "q":
        rows.sort(key=lambda r: (r["q"] if isinstance(r["q"], float)
                                 and math.isfinite(r["q"]) else 9e99))
    elif args.sort == "or":
        rows.sort(key=lambda r: -(r["odds_ratio"]
                                  if math.isfinite(r["odds_ratio"]) else -9e99))

    if args.columns:
        columns = [c.strip() for c in args.columns.split(",")]
    else:
        columns = meta_cols + [c for c in DEFAULT_COLUMNS
                               if c not in meta_cols]

    table_path = os.path.join(args.out, args.table)
    md_path = (os.path.splitext(table_path)[0] + ".md"
               if args.markdown else None)
    print()
    write_table(rows, columns, table_path, md_path)

    full_path = os.path.join(args.out, "titv_full.csv")
    all_fields = list(dict.fromkeys(
        [k for r in rows for k in r if k != "flags"] + ["or_95_ci", "flags"]))
    write_table(rows, all_fields, full_path)

    if args.per_structure:
        for name, mut in muts.items():
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
            with open(os.path.join(args.out, f"{safe}_mutations.csv"),
                      "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(mut[0]),
                                   lineterminator="\n")
                w.writeheader()
                w.writerows(mut)
        print(f"Wrote {len(muts)} per-mutation CSVs to {args.out}/")

    # ---- summary --------------------------------------------------------
    clean = [r for r in rows if not r["flags"]]
    loopy = [r for r in rows if r["direction"] == "loop-biased"]
    print(f"\n{len(clean)} of {len(rows)} structures pass on every "
          f"criterion (q < {args.q_cut}, CI excludes 1, "
          f">= {args.min_subs} substitutions, no empty cell).")
    counts = Counter(f for r in rows for f in r["flags"].split("; ") if f)
    for f, n in counts.most_common():
        print(f"    {f:12s} {n}")
    if loopy:
        print(f"\n  {len(loopy)} structure(s) are loop-biased, the opposite "
              f"of the prediction: "
              f"{', '.join(r['Structure'] for r in loopy[:8])}")
        print(f"  Keep them in the table. Reporting only the rows that "
              f"support the hypothesis is what makes a table look selective.")

    # pooled test across everything analysed
    A = sum(r["stem_ti"] for r in rows)
    B = sum(r["stem_tv"] for r in rows)
    C = sum(r["loop_ti"] for r in rows)
    D = sum(r["loop_tv"] for r in rows)
    if min(A, B, C, D) > 0:
        pooled_or = odds_ratio(A, B, C, D)
        pooled_p = fisher_exact_two_sided(A, B, C, D)
        lo, hi, method = or_ci_bootstrap((A, B, C, D), args.boot, args.seed)
        print(f"\nPOOLED over {len(rows)} structures: "
              f"stem ti/tv {A/B:.2f}, loop ti/tv {C/D:.2f}, "
              f"OR {pooled_or:.2f} ({lo:.2f}-{hi:.2f}), "
              f"Fisher P {fmt_p(pooled_p)}")
        print("  Note: this pools substitutions across related sequences, "
              "which are not independent observations, so P is "
              "anti-conservative. Use --weight and say so in the methods.")
    if problems:
        print(f"\n{len(problems)} input(s) could not be analysed.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
