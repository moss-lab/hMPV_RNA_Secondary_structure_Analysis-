#!/usr/bin/env python3
"""Build the motif summary table from a ScanFold output folder."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

# Filename parsing

# Anchored so motif_5 never matches motif_50 / motif_51.
MOTIF_RE = re.compile(r"motif[_-](\d+)(?=[._\-]|$)", re.IGNORECASE)

# Same, but for free text, capturing the series prefix as well.
TOKEN_RE = re.compile(
    r"([A-Za-z0-9._\-]*?)motif[_-](\d+)(?=[._\-\s]|$)", re.IGNORECASE)

COORD_RE = re.compile(r"coordinates[:_](\d+)[-_](\d+)", re.IGNORECASE)
ZSCORE_RE = re.compile(r"z[-_]?score[=:_](-?\d+(?:\.\d+)?)", re.IGNORECASE)

# Renderings and backups: never a source of numbers.
IGNORED_SUFFIXES = (".bak", ".ps", ".svg", ".pdf", ".png", ".jpg",
                    ".eps", ".log")


def parse_name(basename: str) -> tuple[Optional[str], Optional[int]]:
    """(series, motif number) from a filename, or (None, None)."""
    m = MOTIF_RE.search(basename)
    if not m:
        return None, None
    series = basename[:m.start()].rstrip("_-.") or "(unnamed)"
    return series, int(m.group(1))


def classify(basename: str) -> Optional[str]:
    """The role a file plays. None means ignore it."""
    low = basename.lower()
    if low.endswith(IGNORED_SUFFIXES):
        return None
    if ".r2r." in low or ".dplot." in low or ".surv" in low:
        return None

    if low.endswith(".dbn"):
        return "dbn"
    if low.endswith(".stockholm"):
        return "stockholm"
    if low.endswith(".cm"):
        return "cm"

    # Power files: the covariation source.
    if low.endswith(".cacofold.power"):
        return "power_cacofold"
    if low.endswith(".power"):
        return "power_given"

    if low.endswith(".sorted.cacofold.cov"):
        return "cov_cacofold_sorted"
    if low.endswith(".sorted.cov"):
        return "cov_given_sorted"
    if low.endswith(".cacofold.cov"):
        return "cov_cacofold"
    if low.endswith(".cov"):
        return "cov_given"

    if low.endswith(".cacofold.helixcov"):
        return "helixcov_cacofold"
    if low.endswith(".helixcov"):
        return "helixcov"
    if low.endswith(".cacofold.sto") or low.endswith(".sto"):
        return "sto"
    return None


# Which roles may supply each field, in order of preference.
FIELD_SOURCES = {
    "covariation": {
        "cacofold-power": ["power_cacofold", "power_given"],
        "given-power":    ["power_given", "power_cacofold"],
        "cacofold-cov":   ["cov_cacofold", "cov_cacofold_sorted"],
        "given-cov":      ["cov_given", "cov_given_sorted"],
    },
    "depth": ["stockholm", "power_cacofold", "power_given",
              "cov_given", "cov_cacofold"],
}


# R-scape power / cov parsing

NSEQ_RE = re.compile(r"\bnseq\s+(\d+)", re.IGNORECASE)
ALEN_RE = re.compile(r"\balen\s+(\d+)", re.IGNORECASE)

OBSERVED_RE = re.compile(
    r"BPAIRS\s+observed\s+to\s+covary\s*[:=]?\s*(\d+)", re.IGNORECASE)
EXPECTED_RE = re.compile(
    r"BPAIRS\s+expected\s+to\s+covary\s*[:=]?\s*(\d+(?:\.\d+)?)",
    re.IGNORECASE)
TOTAL_BP_RE = re.compile(r"^#\s*BPAIRS\s+(\d+)\s*$", re.IGNORECASE)

# A data row marked as covarying. R-scape flags these with a leading *,
# and some builds put the marker at the end of the row instead.
STAR_LEAD_RE = re.compile(r"^\s*\*")
STAR_ANY_RE = re.compile(r"(?:^|\s)\*(?:\s|$)")


def parse_power(path: str) -> dict:
    """Covarying base pairs from an R-scape .power file."""
    out = {
        "covarying_bp": None,      # the number used
        "observed_header": None,   # from the header line
        "starred_rows": 0,         # rows carrying a * marker
        "data_rows": 0,            # total data rows (base pairs considered)
        "expected": None,
        "total_bp_header": None,
        "nseq": None,
        "alen": None,
        "basis": None,             # which of the two was used
        "lines": [],               # header lines, for --explain
    }
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                raw = line.rstrip("\n")
                if raw.startswith("#"):
                    if len(out["lines"]) < 40:
                        out["lines"].append(raw)
                    m = OBSERVED_RE.search(raw)
                    if m:
                        out["observed_header"] = int(m.group(1))
                    m = EXPECTED_RE.search(raw)
                    if m:
                        out["expected"] = float(m.group(1))
                    m = TOTAL_BP_RE.match(raw)
                    if m:
                        out["total_bp_header"] = int(m.group(1))
                    if out["nseq"] is None:
                        m = NSEQ_RE.search(raw)
                        if m:
                            out["nseq"] = int(m.group(1))
                    if out["alen"] is None:
                        m = ALEN_RE.search(raw)
                        if m:
                            out["alen"] = int(m.group(1))
                    continue
                if not raw.strip():
                    continue
                out["data_rows"] += 1
                if STAR_LEAD_RE.match(raw) or STAR_ANY_RE.search(raw):
                    out["starred_rows"] += 1
    except OSError as exc:
        print(f"  ! cannot read {path}: {exc}", file=sys.stderr)
        return out

    if out["observed_header"] is not None:
        out["covarying_bp"] = out["observed_header"]
        out["basis"] = "header 'BPAIRS observed to covary'"
    else:
        out["covarying_bp"] = out["starred_rows"]
        out["basis"] = "count of * marked rows (no header line present)"
    return out


def parse_cov(path: str) -> dict:
    """Same quantities from a .cov file, used only as a fallback."""
    return parse_power(path)


def parse_stockholm(path: str) -> dict:
    """Sequence count and alignment length from a Stockholm alignment. Counts."""
    names, alen = set(), None
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("//"):
                    continue
                if line.startswith("#"):
                    continue
                parts = line.split(None, 1)
                if len(parts) == 2:
                    names.add(parts[0])
                    if alen is None:
                        alen = len(parts[1].strip())
    except OSError as exc:
        print(f"  ! cannot read {path}: {exc}", file=sys.stderr)
    return {"nseq": len(names) or None, "alen": alen}


# ScanFold .dbn and summary table

MFE_KEY_RE = re.compile(r"\bMFE\s*[=:]\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE)
ED_KEY_RE = re.compile(
    r"\b(?:ED|ensemble[ _]diversity)\s*[=:]\s*(-?\d+(?:\.\d+)?)",
    re.IGNORECASE)
TRAILING_E_RE = re.compile(r"\(\s*(-?\d+\.\d+)\s*\)\s*$")
SEQ_RE = re.compile(r"^[ACGUTNRYKMSWBDHVacgutnrykmswbdhv]+$")


def parse_dbn(path: str) -> dict:
    out = {"seq": None, "structure": None, "mfe": None, "zscore": None,
           "ed": None, "header": None}
    try:
        with open(path, "r", errors="replace") as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    except OSError as exc:
        print(f"  ! cannot read {path}: {exc}", file=sys.stderr)
        return out
    if not lines:
        return out

    header = lines[0] if lines[0].startswith(">") else ""
    out["header"] = header or None
    body = lines[1:] if header else lines

    for regex, key in ((MFE_KEY_RE, "mfe"), (ED_KEY_RE, "ed"),
                       (ZSCORE_RE, "zscore")):
        m = regex.search(header)
        if m:
            out[key] = float(m.group(1))

    for ln in body:
        token = ln.strip().split(" ")[0]
        if SEQ_RE.match(token):
            if out["seq"] is None:
                out["seq"] = token.upper().replace("T", "U")
        elif token and set(token) <= set("().[]{}<>-_,:"):
            if out["structure"] is None:
                out["structure"] = token
            if out["mfe"] is None:
                m = TRAILING_E_RE.search(ln.strip())
                if m:
                    out["mfe"] = float(m.group(1))
    return out


def parse_scanfold_table(path: str) -> dict:
    """A ScanFold summary table keyed by motif number, with MFE/z/ED."""
    rows: dict[int, dict] = {}
    try:
        with open(path, "r", errors="replace", newline="") as fh:
            sample = fh.read(8192)
            fh.seek(0)
            try:
                reader = csv.reader(fh, csv.Sniffer().sniff(
                    sample, delimiters="\t,;"))
            except csv.Error:
                reader = (re.split(r"\s+", ln.strip()) for ln in fh
                          if ln.strip())
            header = None
            for parts in reader:
                if not parts:
                    continue
                if header is None:
                    header = [p.strip().lower() for p in parts]
                    continue
                rec = dict(zip(header, [p.strip() for p in parts]))
                num = None
                for key, val in rec.items():
                    if "motif" in key and val:
                        m = re.search(r"(\d+)", val)
                        if m:
                            num = int(m.group(1))
                            break
                if num is None:
                    continue

                def pick(*needles):
                    for k, v in rec.items():
                        if any(n in k for n in needles) and v:
                            try:
                                return float(v)
                            except ValueError:
                                return None
                    return None

                rows[num] = {"mfe": pick("mfe", "energy"),
                             "zscore": pick("zscore", "z-score", "z_score"),
                             "ed": pick("ed", "ensemble")}
    except OSError as exc:
        print(f"! cannot read ScanFold table {path}: {exc}", file=sys.stderr)
    return rows


def compute_mfe_ed(seq: str) -> tuple[Optional[float], Optional[float]]:
    """Refold with ViennaRNA if importable."""
    try:
        import RNA  # type: ignore
    except ImportError:
        return None, None
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    return round(mfe, 2), round(fc.mean_bp_distance(), 2)


# All.CaCo.pdf

def extract_pdf_pages(path: str) -> list[str]:
    try:
        import pdfplumber  # type: ignore
        with pdfplumber.open(path) as pdf:
            return [(pg.extract_text() or "") for pg in pdf.pages]
    except ImportError:
        pass
    except Exception as exc:
        print(f"! pdfplumber failed: {exc}", file=sys.stderr)

    PdfReader = None
    try:
        from pypdf import PdfReader  # type: ignore
    except ImportError:
        try:
            from PyPDF2 import PdfReader  # type: ignore
        except ImportError:
            pass
    if PdfReader is not None:
        try:
            return [(p.extract_text() or "") for p in PdfReader(path).pages]
        except Exception as exc:
            print(f"! pypdf failed: {exc}", file=sys.stderr)

    import shutil
    import subprocess
    if shutil.which("pdftotext"):
        try:
            return subprocess.run(["pdftotext", "-layout", path, "-"],
                                  capture_output=True, text=True,
                                  check=True).stdout.split("\f")
        except subprocess.CalledProcessError as exc:
            print(f"! pdftotext failed: {exc}", file=sys.stderr)

    sys.exit(
        "Cannot read the PDF: no text extractor available.\n"
        "    pip install --user pypdf\n"
        "or extract the text first and pass it in:\n"
        "    pdftotext -layout All.CaCo.pdf caco.txt\n"
        "    python3 build_motif_table.py --root ... --pdf-text caco.txt")


def read_pdf_text_file(path: str) -> list[str]:
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    if "\f" in text:
        return text.split("\f")
    if "===== page " in text:
        return [p for p in re.split(r"=====\s*page\s+\d+\s*=====", text)
                if p.strip()]
    return [text]


def parse_caco_pdf(path: str, dump: Optional[str] = None,
                   text_file: Optional[str] = None) -> dict[tuple, dict]:
    """Membership and coordinates per motif. The pages carry a header line and."""
    pages = read_pdf_text_file(text_file) if text_file \
        else extract_pdf_pages(path)
    if dump:
        with open(dump, "w") as fh:
            for i, text in enumerate(pages, 1):
                fh.write(f"===== page {i} =====\n{text}\n")
        print(f"Wrote extracted PDF text to {dump}")

    found: dict[tuple, dict] = {}
    for i, text in enumerate(pages, 1):
        if not text.strip():
            continue
        keys = {((m.group(1) or "").rstrip("_-.") or "(unnamed)",
                 int(m.group(2))) for m in TOKEN_RE.finditer(text)}
        if not keys:
            continue
        if len(keys) > 1:
            print(f"  ! PDF page {i} names several motifs {sorted(keys)}; "
                  f"skipped", file=sys.stderr)
            continue
        key = keys.pop()
        rec = found.setdefault(key, {"pages": [], "kinds": set()})
        rec["pages"].append(i)
        head = text.split("\n")[0]
        rec["kinds"].add("cacofold" if "cacofold" in head.lower() else "given")
        if rec.get("zscore") is None:
            m = ZSCORE_RE.search(text)
            if m:
                rec["zscore"] = float(m.group(1))
        if rec.get("start") is None:
            m = COORD_RE.search(text)
            if m:
                rec["start"], rec["end"] = int(m.group(1)), int(m.group(2))
    return found


# Assembly

@dataclass
class Motif:
    number: int
    series: str = ""
    files: dict = field(default_factory=lambda: defaultdict(list))
    start: Optional[int] = None
    end: Optional[int] = None
    zscore_fn: Optional[float] = None

    def best(self, role: str) -> Optional[str]:
        paths = self.files.get(role) or []
        return max(paths, key=os.path.getmtime) if paths else None

    def first_of(self, roles) -> tuple[Optional[str], Optional[str]]:
        for role in roles:
            path = self.best(role)
            if path:
                return role, path
        return None, None


def collect(root: str) -> dict[tuple, Motif]:
    motifs: dict[tuple, Motif] = {}
    for dirpath, _dirs, filenames in os.walk(root):
        for name in filenames:
            series, num = parse_name(name)
            if num is None:
                continue
            role = classify(name)
            if role is None:
                continue
            key = (series, num)
            mot = motifs.setdefault(key, Motif(number=num, series=series))
            mot.files[role].append(os.path.join(dirpath, name))

            c = COORD_RE.search(name)
            if c:
                start, end = int(c.group(1)), int(c.group(2))
                if mot.start is None:
                    mot.start, mot.end = start, end
                elif (mot.start, mot.end) != (start, end):
                    print(f"  ! {series} motif {num}: files disagree on "
                          f"coordinates, {mot.start}-{mot.end} vs "
                          f"{start}-{end} ({name})", file=sys.stderr)
            z = ZSCORE_RE.search(name)
            if z and mot.zscore_fn is None:
                mot.zscore_fn = float(z.group(1))
    return motifs


def resolve_pdf_keys(pdf_rows, motifs, args):
    """Tie each PDF entry to a motif on disk. Coordinates are the only."""
    on_disk = set(motifs)
    by_coord = defaultdict(list)
    for key, mot in motifs.items():
        if mot.start is not None:
            by_coord[(mot.start, mot.end)].append(key)

    resolved, n_coord, n_num = {}, 0, 0
    for (p_series, p_num), rec in pdf_rows.items():
        if args.match_by in ("auto", "number") and (p_series, p_num) in on_disk:
            resolved[(p_series, p_num)] = rec
            continue

        span = (rec.get("start"), rec.get("end"))
        hits = by_coord.get(span, []) if span[0] is not None else []
        if args.match_by in ("auto", "coords") and len(hits) == 1:
            if hits[0][1] != p_num:
                print(f"  ! PDF motif {p_num} ({p_series}) matched to motif "
                      f"{hits[0][1]} ({hits[0][0]}) via coordinates "
                      f"{span[0]}-{span[1]}", file=sys.stderr)
            resolved[hits[0]] = rec
            n_coord += 1
            continue
        if len(hits) > 1:
            print(f"  ! coordinates {span[0]}-{span[1]} match {len(hits)} "
                  f"motifs; skipped", file=sys.stderr)
            continue
        if args.match_by == "coords":
            print(f"  ! PDF motif {p_num} ({p_series}): nothing on disk at "
                  f"{span[0]}-{span[1]}; skipped", file=sys.stderr)
            continue

        # Number matching is only safe inside one series.
        same = [k for k in on_disk if k[1] == p_num
                and (k[0] == p_series or p_series == "(unnamed)")]
        cross = [k for k in on_disk if k[1] == p_num and k not in same]
        if not same and cross:
            print(f"  ! PDF motif {p_num} ({p_series}) has no coordinate "
                  f"match, and motif {p_num} of {cross[0][0]!r} is a "
                  f"different structure. Skipped rather than guessed.",
                  file=sys.stderr)
        if len(same) == 1:
            resolved[same[0]] = rec
            n_num += 1
        elif len(same) > 1:
            print(f"  ! PDF motif {p_num} matches {len(same)} motifs on "
                  f"disk; use --series", file=sys.stderr)
        else:
            resolved[(p_series, p_num)] = rec
    if n_coord or n_num:
        print(f"Matched PDF entries: {n_coord} by coordinates, "
              f"{n_num} by motif number")
    return resolved


def build_rows(motifs, args, scanfold_rows, pdf_rows):
    rows, audit, warnings = [], [], []
    multi_series = len({k[0] for k in motifs}) > 1
    cov_roles = FIELD_SOURCES["covariation"][args.covariation]

    for key in sorted(motifs, key=lambda k: (k[0], k[1])):
        mot, num = motifs[key], motifs[key].number
        src = {"series": mot.series}
        pdf = pdf_rows.get(key, {})

        if pdf:
            src["pdf_pages"] = ",".join(str(p) for p in pdf.get("pages", []))
            if mot.start is None and pdf.get("start") is not None:
                mot.start, mot.end = pdf["start"], pdf["end"]
            elif (pdf.get("start") is not None
                  and (mot.start, mot.end) != (pdf["start"], pdf["end"])):
                warnings.append(
                    f"motif {num}: coordinates {mot.start}-{mot.end} on disk "
                    f"vs {pdf['start']}-{pdf['end']} in the PDF")

        # ---- structure metrics ------------------------------------------
        dbn_path = mot.best("dbn")
        dbn = parse_dbn(dbn_path) if dbn_path else {}
        if dbn_path:
            src["dbn"] = dbn_path

        mfe, ed = dbn.get("mfe"), dbn.get("ed")
        zscore = dbn.get("zscore")
        if mfe is not None:
            src["mfe"] = dbn_path
        if ed is not None:
            src["ed"] = dbn_path

        rec = scanfold_rows.get(num, {})
        if mfe is None and rec.get("mfe") is not None:
            mfe, src["mfe"] = rec["mfe"], args.scanfold_table
        if ed is None and rec.get("ed") is not None:
            ed, src["ed"] = rec["ed"], args.scanfold_table
        if zscore is None:
            zscore = rec.get("zscore")

        if zscore is None and pdf.get("zscore") is not None:
            zscore, src["zscore"] = pdf["zscore"], "All.CaCo.pdf"
        if zscore is None:
            zscore = mot.zscore_fn
            if zscore is not None:
                src["zscore"] = "filename"
        elif mot.zscore_fn is not None and abs(zscore - mot.zscore_fn) > 0.011:
            warnings.append(f"motif {num}: z-score {zscore} in file contents "
                            f"vs {mot.zscore_fn} in the filename")

        if (mfe is None or ed is None) and args.compute_missing and dbn.get("seq"):
            c_mfe, c_ed = compute_mfe_ed(dbn["seq"])
            if mfe is None and c_mfe is not None:
                mfe, src["mfe"] = c_mfe, "recomputed (ViennaRNA)"
            if ed is None and c_ed is not None:
                ed, src["ed"] = c_ed, "recomputed (ViennaRNA)"

        # ---- covarying base pairs ---------------------------------------
        cov_bp, power = None, {}
        role, cov_path = mot.first_of(cov_roles)
        if cov_path:
            power = parse_power(cov_path)
            cov_bp = power["covarying_bp"]
            src["covarying_bp"] = cov_path
            src["covarying_basis"] = power["basis"]
            if power["observed_header"] is None:
                warnings.append(
                    f"motif {num}: no 'BPAIRS observed to covary' line in "
                    f"{os.path.basename(cov_path)}; counted "
                    f"{power['starred_rows']} marked rows instead")
            elif (power["starred_rows"]
                  and power["starred_rows"] != power["observed_header"]):
                warnings.append(
                    f"motif {num}: header says {power['observed_header']} "
                    f"covarying pairs but {power['starred_rows']} rows are "
                    f"marked ({os.path.basename(cov_path)})")
        else:
            warnings.append(f"motif {num}: no {args.covariation} file found")

        # ---- alignment depth --------------------------------------------
        depth = None
        for role in FIELD_SOURCES["depth"]:
            path = mot.best(role)
            if not path:
                continue
            parsed = (parse_stockholm(path) if role == "stockholm"
                      else parse_power(path))
            if parsed.get("nseq"):
                depth = parsed["nseq"]
                src["depth"] = path
                break

        # Sanity: significant covariation needs sequences to covary in.
        if cov_bp and depth is not None and cov_bp > depth:
            warnings.append(
                f"motif {num}: {cov_bp} covarying pairs reported from an "
                f"alignment of {depth} sequences, which is not possible. "
                f"Check {os.path.basename(src.get('covarying_bp', '?'))}")

        # ---- length ------------------------------------------------------
        length = (mot.end - mot.start + 1) if mot.start is not None else None
        if length is None and dbn.get("seq"):
            length = len(dbn["seq"])

        row = {}
        if multi_series:
            row["Series"] = mot.series
        row.update({
            "Motif no.": num,
            "Coordinates": f"{mot.start}-{mot.end}" if mot.start else "",
            "MFE (kcal/mol)": "" if mfe is None else mfe,
            "z-score": "" if zscore is None else zscore,
            "ED": "" if ed is None else ed,
            "Covarying bp": "" if cov_bp is None else cov_bp,
            "Alignment depth (n)": "" if depth is None else depth,
            "Length_Motif": "" if length is None else length,
        })
        if args.extra_columns:
            row["BP considered"] = power.get("data_rows", "")
            row["Expected to covary"] = (power.get("expected")
                                         if power.get("expected") is not None
                                         else "")
        rows.append(row)
        audit.append({"Motif no.": num, **src})

    if rows:
        cols = list(rows[0].keys())
        for r in rows:
            for c in cols:
                r.setdefault(c, "")
    return rows, audit, warnings


# Output

# Decimal places per column. Anything not listed is written as-is, and
# integers stay integers rather than gaining a spurious ".0".
DECIMALS = {
    "MFE (kcal/mol)": 2,
    "z-score": 2,
    "ED": 2,
    "Expected to covary": 2,
}

# Fixed column order, so both tables and any rerun agree.
COLUMN_ORDER = [
    "Series", "Motif no.", "Coordinates", "MFE (kcal/mol)", "z-score", "ED",
    "Covarying bp", "Alignment depth (n)", "Length_Motif",
    "BP considered", "Expected to covary",
]


def ordered_columns(rows):
    present = set(rows[0].keys())
    cols = [c for c in COLUMN_ORDER if c in present]
    cols += [c for c in rows[0].keys() if c not in COLUMN_ORDER]
    return cols


def format_cell(col, value, na="NA"):
    """One rendering rule per column, so the file is consistent throughout."""
    if value is None or value == "":
        return na
    places = DECIMALS.get(col)
    if places is not None:
        try:
            return f"{float(value):.{places}f}"
        except (TypeError, ValueError):
            return str(value)
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def write_table(rows, path, markdown=None, na="NA", label=""):
    if not rows:
        print("No motifs to write.", file=sys.stderr)
        return
    cols = ordered_columns(rows)
    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=delim, lineterminator="\n",
                       quoting=csv.QUOTE_MINIMAL)
        w.writerow(cols)
        for r in rows:
            w.writerow([format_cell(c, r.get(c, ""), na) for c in cols])
    what = f"{label} " if label else ""
    print(f"Wrote {len(rows)} motifs to {what}{path}")

    if markdown:
        widths = [max(len(c), *(len(format_cell(c, r.get(c, ""), na))
                                for r in rows)) for c in cols]
        numeric = {c for c in cols if c in DECIMALS
                   or c in ("Motif no.", "Covarying bp",
                            "Alignment depth (n)", "Length_Motif",
                            "BP considered")}
        def cell(c, text, w):
            return text.rjust(w) if c in numeric else text.ljust(w)
        with open(markdown, "w") as fh:
            fh.write("| " + " | ".join(cell(c, c, w)
                                       for c, w in zip(cols, widths)) + " |\n")
            fh.write("|" + "|".join(("-" * (w + 2)) if c not in numeric
                                    else ("-" * (w + 1) + ":")
                                    for c, w in zip(cols, widths)) + "|\n")
            for r in rows:
                fh.write("| " + " | ".join(
                    cell(c, format_cell(c, r.get(c, ""), na), w)
                    for c, w in zip(cols, widths)) + " |\n")
        print(f"Wrote markdown to {markdown}")


AUDIT_ORDER = ["Motif no.", "series", "covarying_bp", "covarying_basis",
               "depth", "mfe", "ed", "zscore", "dbn", "pdf_pages"]


def write_audit(audit, path, na="NA", root=None):
    """One row per motif naming the file behind each value. Paths are written."""
    keys = sorted({k for row in audit for k in row})
    cols = [c for c in AUDIT_ORDER if c in keys]
    cols += [c for c in keys if c not in AUDIT_ORDER]

    def shorten(v):
        if root and isinstance(v, str) and v.startswith(root):
            return os.path.relpath(v, root)
        return v

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for row in sorted(audit, key=lambda r: r["Motif no."]):
            w.writerow([shorten(row.get(c, "")) if row.get(c, "") != ""
                        else na for c in cols])
    print(f"Wrote per-field sources to {path}")


def write_crosswalk(all_motifs, path, na="NA"):
    series = sorted({k[0] for k in all_motifs})
    by_span = {}
    for key, mot in all_motifs.items():
        if mot.start is None:
            continue
        by_span.setdefault((mot.start, mot.end), {})[key[0]] = mot.number
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Start", "End", "Length"] + series)
        for span in sorted(by_span):
            row = by_span[span]
            w.writerow([span[0], span[1], span[1] - span[0] + 1]
                       + [row.get(s, na) for s in series])
    shared = sum(1 for v in by_span.values() if len(v) == len(series))
    print(f"Wrote crosswalk to {path} ({len(by_span)} spans, "
          f"{shared} in every series)")


def explain(motifs, num, args):
    """Show how one motif's covariation number was derived."""
    targets = [k for k in motifs if k[1] == num]
    if not targets:
        print(f"No motif {num} found.")
        return
    for key in targets:
        mot = motifs[key]
        print(f"\n=== [{key[0]}] motif {num}  "
              f"{mot.start}-{mot.end} ===")
        for role in ("power_cacofold", "power_given", "cov_cacofold",
                     "cov_given", "stockholm", "dbn"):
            path = mot.best(role)
            if path:
                print(f"  {role:18s} {os.path.basename(path)}")
        role, path = mot.first_of(FIELD_SOURCES["covariation"][args.covariation])
        if not path:
            print(f"  no {args.covariation} file present")
            continue
        p = parse_power(path)
        print(f"\n  Using: {os.path.basename(path)}")
        print(f"  header 'observed to covary' : {p['observed_header']}")
        print(f"  rows marked with *          : {p['starred_rows']}")
        print(f"  total data rows             : {p['data_rows']}")
        print(f"  expected to covary          : {p['expected']}")
        print(f"  nseq / alen in header       : {p['nseq']} / {p['alen']}")
        print(f"  -> Covarying bp = {p['covarying_bp']}  ({p['basis']})")
        sto = mot.best("stockholm")
        if sto:
            s = parse_stockholm(sto)
            print(f"  sequences in .stockholm     : {s['nseq']}")
        print("\n  --- header lines of that file ---")
        for ln in p["lines"][:25]:
            print("   ", ln)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root",
                    default="/lustre/hdd/LAS/wmoss-lab/ddabiri/CaCo_folder/MNV")
    ap.add_argument("--out", default="motif_table.tsv")
    ap.add_argument("--markdown")
    ap.add_argument("--audit")
    ap.add_argument("--crosswalk")
    ap.add_argument("--series", action="append",
                    help="keep only this ScanFold run, e.g. "
                         "--series NC_hMNV_cRNA. Repeatable.")
    ap.add_argument("--covariation", default="cacofold-power",
                    choices=["cacofold-power", "given-power",
                             "cacofold-cov", "given-cov"],
                    help="where Covarying bp comes from "
                         "(default: cacofold-power)")
    ap.add_argument("--caco-pdf",
                    help="defaults to <root>/All.CaCo.pdf when present")
    ap.add_argument("--pdf-text",
                    help="use pre-extracted PDF text instead of the PDF")
    ap.add_argument("--dump-pdf-text")
    ap.add_argument("--no-restrict", action="store_true",
                    help="keep every motif on disk, not only those in the PDF")
    ap.add_argument("--match-by", choices=["auto", "coords", "number"],
                    default="auto")
    ap.add_argument("--scanfold-table",
                    help="ScanFold summary table carrying MFE/z-score/ED")
    ap.add_argument("--compute-missing", action="store_true",
                    help="refold with ViennaRNA to fill missing MFE/ED")
    ap.add_argument("--covarying-out", nargs="?", const="AUTO",
                    metavar="FILE",
                    help="also write a table of only the motifs with "
                         "covarying base pairs. With no value, the name is "
                         "derived from --out (motif_table_covarying.tsv).")
    ap.add_argument("--min-covarying", type=int, default=1,
                    help="how many covarying pairs a motif needs to appear "
                         "in that file (default 1)")
    ap.add_argument("--na", default="NA",
                    help="what to write for a missing value "
                         "(default NA; pass '' for an empty cell)")
    ap.add_argument("--extra-columns", action="store_true",
                    help="also emit BP considered and Expected to covary")
    ap.add_argument("--explain", type=int, metavar="MOTIF",
                    help="show how this motif's covariation was derived, "
                         "then exit")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    if not os.path.isdir(args.root):
        sys.exit(f"Not a directory: {args.root}")

    motifs = collect(args.root)
    all_motifs = dict(motifs)

    counts = defaultdict(int)
    for (series, _n) in motifs:
        counts[series] += 1
    print(f"Found {len(motifs)} motif records under {args.root}")
    for series, n in sorted(counts.items()):
        print(f"    series {series!r}: {n} motifs")

    if args.series:
        keep = set(args.series)
        unknown = keep - set(counts)
        if unknown:
            sys.exit(f"No such series: {', '.join(sorted(unknown))}. "
                     f"Available: {', '.join(sorted(counts))}")
        motifs = {k: v for k, v in motifs.items() if k[0] in keep}
        print(f"Restricted to {', '.join(sorted(keep))}: {len(motifs)} motifs")
    elif len(counts) > 1:
        print("\n!! More than one ScanFold run is present and each numbers "
              "its motifs independently. Use --series to pick the one your "
              "manuscript uses.\n", file=sys.stderr)

    if args.explain is not None:
        explain(motifs, args.explain, args)
        return

    if args.crosswalk:
        write_crosswalk(all_motifs, args.crosswalk, args.na)

    # ---- PDF membership -------------------------------------------------
    pdf_rows = {}
    pdf_source = args.pdf_text or args.caco_pdf
    if pdf_source is None:
        default_pdf = os.path.join(args.root, "All.CaCo.pdf")
        if os.path.isfile(default_pdf):
            args.caco_pdf, pdf_source = default_pdf, default_pdf
    if pdf_source:
        if not os.path.isfile(pdf_source):
            sys.exit(f"No such file: {pdf_source}")
        pdf_rows = parse_caco_pdf(args.caco_pdf or pdf_source,
                                  args.dump_pdf_text, args.pdf_text)
        pdf_rows = resolve_pdf_keys(pdf_rows, motifs, args)

        in_pdf, on_disk = set(pdf_rows), set(motifs)
        print(f"{os.path.basename(pdf_source)} lists {len(in_pdf)} motifs")

        only_given = sorted(k[1] for k, v in pdf_rows.items()
                            if v.get("kinds") == {"given"})
        if only_given:
            print(f"  note: no CaCoFold page for motifs "
                  f"{', '.join(map(str, only_given))}")

        missing = sorted(in_pdf - on_disk, key=lambda k: k[1])
        kept_series = set(args.series) if args.series else None
        for key in missing:
            if kept_series is not None and key[0] not in kept_series:
                print(f"  ! PDF motif {key[1]} is in series {key[0]!r}, "
                      f"excluded by --series", file=sys.stderr)
                continue
            print(f"  ! PDF motif {key[1]}: no files on disk",
                  file=sys.stderr)
            motifs.setdefault(key, Motif(number=key[1], series=key[0]))

        if not args.no_restrict:
            dropped = sorted(on_disk - in_pdf, key=lambda k: k[1])
            motifs = {k: v for k, v in motifs.items() if k in in_pdf}
            if dropped:
                print(f"Excluded {len(dropped)} motifs absent from the PDF: "
                      f"{', '.join(str(d[1]) for d in dropped)}")

    if args.debug:
        for key in sorted(motifs, key=lambda k: (k[0], k[1])):
            print(f"\n[{key[0]}] motif {key[1]}:")
            for role, paths in sorted(motifs[key].files.items()):
                for p in paths:
                    print(f"   {role:20s} {os.path.basename(p)}")

    scanfold_rows = (parse_scanfold_table(args.scanfold_table)
                     if args.scanfold_table else {})

    rows, audit, warnings = build_rows(motifs, args, scanfold_rows, pdf_rows)
    write_table(rows, args.out, args.markdown, args.na,
                label="all motifs ->")
    if args.audit:
        write_audit(audit, args.audit, args.na, args.root)

    if args.covarying_out:
        path = args.covarying_out
        if path == "AUTO":
            stem, ext = os.path.splitext(args.out)
            path = f"{stem}_covarying{ext or '.tsv'}"
        cov_rows, no_data = [], []
        for r in rows:
            val = r["Covarying bp"]
            if val == "":
                no_data.append(r["Motif no."])
            elif int(val) >= args.min_covarying:
                cov_rows.append(r)
        md = None
        if args.markdown:
            stem, ext = os.path.splitext(args.markdown)
            md = f"{stem}_covarying{ext or '.md'}"
        if cov_rows:
            write_table(cov_rows, path, md, args.na,
                        label="covarying motifs ->")
            total = sum(int(r["Covarying bp"]) for r in cov_rows)
            print(f"    {len(cov_rows)} of {len(rows)} motifs carry at least "
                  f"{args.min_covarying} covarying bp "
                  f"({total} covarying pairs in total)")
        else:
            print(f"No motif reached {args.min_covarying} covarying bp; "
                  f"{path} not written", file=sys.stderr)
        if no_data:
            print(f"    ! {len(no_data)} motifs had no covariation file and "
                  f"are in neither group: "
                  f"{', '.join(str(n) for n in no_data)}", file=sys.stderr)

    if warnings:
        print(f"\n{len(warnings)} warnings:", file=sys.stderr)
        for w in warnings:
            print(f"  ! {w}", file=sys.stderr)

    blank_mfe = [r["Motif no."] for r in rows if r["MFE (kcal/mol)"] == ""]
    blank_ed = [r["Motif no."] for r in rows if r["ED"] == ""]
    if blank_mfe or blank_ed:
        print(f"\nMFE missing for {len(blank_mfe)} motifs, ED missing for "
              f"{len(blank_ed)}. Neither is in All.CaCo.pdf. Point "
              f"--scanfold-table at the ScanFold output that carries them, "
              f"or use --compute-missing to refold.")


if __name__ == "__main__":
    main()
