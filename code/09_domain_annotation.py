"""
09_domain_annotation.py
=======================
UniProt feature mining for the Dimer-Gain (Interface Unmasking) cohort.


Run order: 9 of 14
"""

import re
import time
from pathlib import Path

import pandas as pd
import requests

try:
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment  # noqa: F401
except ImportError:
    import subprocess, sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "biopython"])
    from Bio import pairwise2


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
INPUT_CSV     = "Dimer_Gain_Cohort.csv"
OUT_LONG      = "Domain_Annotations_v2.csv"
OUT_PER_EVENT = "Domain_Annotations_v2_per_event.csv"
OUT_SUMMARY   = "Domain_Annotations_v2_summary.txt"

UNIPROT_JSON  = "https://rest.uniprot.org/uniprotkb/{acc}.json"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

FEATURE_TYPES_KEEP = {
    "Domain", "Region", "Motif", "Repeat", "Compositional bias",
    "Coiled coil", "Zinc finger", "DNA binding", "Active site",
    "Topological domain", "Transmembrane", "Signal peptide", "Propeptide",
    "Binding site",
}

INHIBITORY_PATTERNS = re.compile(
    r"\b(autoinhibit|auto-inhibit|inhibitory|pseudo[- ]?substrate|"
    r"regulatory|pre[- ]?segment|prodomain|propeptide|"
    r"intramolecular[- ]?lock|masking|occlud)",
    re.IGNORECASE,
)
DISORDER_PATTERNS = re.compile(
    r"\b(disordered|low complexity|polar|acidic|basic|gly[- ]rich|pro[- ]rich|"
    r"ser[- ]rich|polyampholyte|compositional bias)",
    re.IGNORECASE,
)

REQUEST_HEADERS = {
    "User-Agent": "KaragolEtAl-Manuscript/1.0",
    "Accept":     "application/json",
}


# ---------------------------------------------------------------------------
# UNIPROT FETCH HELPERS
# ---------------------------------------------------------------------------
def _get(url, accept="application/json", retries=3, sleep=0.2):
    """Polite GET with exponential backoff."""
    for k in range(retries):
        try:
            r = requests.get(url, headers={**REQUEST_HEADERS, "Accept": accept},
                             timeout=30)
            if r.status_code == 200:
                return r.text if accept != "application/json" else r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 ** k)
                continue
            return None
        except requests.RequestException:
            time.sleep(2 ** k)
    return None


def fetch_canonical_entry(acc):
    return _get(UNIPROT_JSON.format(acc=acc), accept="application/json")


def fetch_isoform_fasta(iso_acc):
    txt = _get(UNIPROT_FASTA.format(acc=iso_acc), accept="text/x-fasta")
    if not txt:
        return None
    return "".join(line.strip() for line in txt.splitlines()
                   if not line.startswith(">")) or None


def extract_features(entry):
    """Return list of dicts: {type, description, begin, end}."""
    out = []
    for f in entry.get("features", []):
        ftype = f.get("type")
        if ftype not in FEATURE_TYPES_KEEP:
            continue
        loc = f.get("location", {}) or {}
        try:
            begin = int(loc.get("start", {}).get("value"))
            end   = int(loc.get("end",   {}).get("value"))
        except (TypeError, ValueError):
            continue
        out.append({
            "type":        ftype,
            "description": f.get("description", "") or "",
            "begin":       begin,
            "end":         end,
        })
    return out


def extract_canonical_seq(entry):
    return entry.get("sequence", {}).get("value") or ""


# ---------------------------------------------------------------------------
# DELETION WINDOW DERIVATION (alignment-based, topology-agnostic)
# ---------------------------------------------------------------------------
def derive_deletion_windows(canon_seq, iso_seq):
    """
    Globally align iso onto canon with affine gap penalties tuned for
    splice-style indels (single large gap rather than many small ones).
    Return list of (start, end) intervals in canonical numbering (1-based,
    inclusive) that are deleted in the isoform, plus a topology label.
    """
    if not canon_seq or not iso_seq:
        return [], "Unknown"

    alns = pairwise2.align.globalms(
        canon_seq, iso_seq, 2, -1, -10, -0.5, one_alignment_only=True
    )
    if not alns:
        return [], "Unknown"
    a_canon, a_iso, _, _, _ = alns[0]

    deletions = []
    canon_pos = 0
    cur_start = None
    cur_end   = None
    for c_char, i_char in zip(a_canon, a_iso):
        if c_char != "-":
            canon_pos += 1
        if c_char != "-" and i_char == "-":
            if cur_start is None:
                cur_start = canon_pos
            cur_end = canon_pos
        else:
            if cur_start is not None:
                deletions.append((cur_start, cur_end))
                cur_start = None
    if cur_start is not None:
        deletions.append((cur_start, cur_end))

    if not deletions:
        return [], "No-Deletion"

    largest = max(deletions, key=lambda x: x[1] - x[0] + 1)
    L = len(canon_seq)
    if largest[0] <= 5:
        topo = "N-Terminal Loss"
    elif largest[1] >= L - 5:
        topo = "C-Terminal Loss"
    else:
        topo = "Internal Deletion / Splice"

    return deletions, topo


def overlap(a, b):
    """Inclusive interval overlap length."""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


# ---------------------------------------------------------------------------
# MAIN PIPELINE
# ---------------------------------------------------------------------------
def main():
    df = pd.read_csv(INPUT_CSV)
    n_events = len(df)
    print(f"Loaded {n_events} (canonical, isoform) events spanning "
          f"{df['gene'].nunique()} unique genes.")

    canon_cache    = {}
    long_rows      = []
    per_event_rows = []

    for idx, row in df.iterrows():
        gene      = row["gene"]
        canon_acc = row["canonical_acc"]
        iso_acc   = row["isoform_acc"]

        if canon_acc not in canon_cache:
            canon_cache[canon_acc] = fetch_canonical_entry(canon_acc)
            time.sleep(0.2)
        entry = canon_cache.get(canon_acc)
        if entry is None:
            print(f"  [WARN] {gene} {canon_acc} - UniProt fetch failed")
            continue

        canon_seq = extract_canonical_seq(entry)
        features  = extract_features(entry)

        iso_seq = fetch_isoform_fasta(iso_acc)
        time.sleep(0.2)
        if iso_seq is None:
            iso_len   = int(row["iso_len"])
            canon_len = int(row["canon_len"])
            deletions = ([(iso_len + 1, canon_len)]
                         if canon_len > iso_len else [])
            topology  = "Unresolved (length-fallback)"
        else:
            deletions, topology = derive_deletion_windows(canon_seq, iso_seq)

        deletion_total_len = sum(e - s + 1 for s, e in deletions)
        deletion_str       = ";".join(f"{s}-{e}" for s, e in deletions) or ""

        overlapping     = []
        any_inhibitory  = False
        any_disorder    = False
        for feat in features:
            ov = sum(overlap((feat["begin"], feat["end"]), d) for d in deletions)
            if ov == 0:
                continue
            feat_len   = feat["end"] - feat["begin"] + 1
            kind       = "full" if ov >= feat_len else "partial"
            inhibitory = bool(INHIBITORY_PATTERNS.search(
                f'{feat["type"]} {feat["description"]}'))
            disorder   = bool(DISORDER_PATTERNS.search(
                f'{feat["type"]} {feat["description"]}'))
            any_inhibitory |= inhibitory
            any_disorder   |= disorder
            long_rows.append({
                "gene":                  gene,
                "canonical_acc":         canon_acc,
                "isoform_acc":           iso_acc,
                "topology":              topology,
                "deletion_window":       deletion_str,
                "deletion_length":       deletion_total_len,
                "feature_type":          feat["type"],
                "feature_description":   feat["description"],
                "feature_begin":         feat["begin"],
                "feature_end":           feat["end"],
                "overlap_residues":      ov,
                "overlap_kind":          kind,
                "is_inhibitory_lexical": inhibitory,
                "is_disorder_lexical":   disorder,
            })
            overlapping.append(
                f'{feat["type"]}:{feat["description"]}'
                f'[{feat["begin"]}-{feat["end"]}]({kind})'
            )

        per_event_rows.append({
            "gene":                    gene,
            "canonical_acc":           canon_acc,
            "isoform_acc":             iso_acc,
            "canon_len":               int(row["canon_len"]),
            "iso_len":                 int(row["iso_len"]),
            "topology":                topology,
            "deletion_window":         deletion_str,
            "deletion_length":         deletion_total_len,
            "n_features_overlapping":  len(overlapping),
            "features_overlapping":    "; ".join(overlapping),
            "any_inhibitory_lexical":  any_inhibitory,
            "any_disorder_lexical":    any_disorder,
        })

        if (idx + 1) % 25 == 0:
            print(f"  ...{idx+1}/{n_events} events processed")

    long_df = pd.DataFrame(long_rows)
    pe_df   = pd.DataFrame(per_event_rows)
    long_df.to_csv(OUT_LONG, index=False)
    pe_df.to_csv(OUT_PER_EVENT, index=False)

    # ---- Summary
    n           = len(pe_df)
    deletes_any = (pe_df["n_features_overlapping"] > 0).sum()
    deletes_inh = pe_df["any_inhibitory_lexical"].sum()
    deletes_dis = pe_df["any_disorder_lexical"].sum()

    type_counts = (long_df["feature_type"].value_counts()
                   if not long_df.empty else pd.Series(dtype=int))
    topo_counts = pe_df["topology"].value_counts()

    summary = []
    summary.append("PATH A - UNIPROT FEATURE MINING SUMMARY")
    summary.append("=" * 60)
    summary.append(f"Events analysed: {n}")
    summary.append(f"Events deleting >=1 curated feature: "
                   f"{deletes_any} ({100*deletes_any/n:.1f}%)")
    summary.append(f"Events deleting an INHIBITORY-annotated feature "
                   f"(lexical): {deletes_inh} ({100*deletes_inh/n:.1f}%)")
    summary.append(f"Events deleting a DISORDER/COMPOSITIONAL-BIAS feature: "
                   f"{deletes_dis} ({100*deletes_dis/n:.1f}%)")
    summary.append("")
    summary.append("Topology distribution (alignment-derived):")
    for k, v in topo_counts.items():
        summary.append(f"  {k}: {v} ({100*v/n:.1f}%)")
    summary.append("")
    summary.append("Most-frequent feature types deleted (top 15):")
    for k, v in type_counts.head(15).items():
        summary.append(f"  {k}: {v}")
    summary.append("")
    summary.append("Outputs:")
    summary.append(f"  {OUT_LONG}        (one row per event x feature)")
    summary.append(f"  {OUT_PER_EVENT}  (one row per event)")

    text = "\n".join(summary)
    Path(OUT_SUMMARY).write_text(text)
    print()
    print(text)


if __name__ == "__main__":
    main()
