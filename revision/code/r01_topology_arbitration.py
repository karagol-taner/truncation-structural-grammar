"""Arbitrate two disagreeing topology assignments against an external criterion.

Two of the cohort tables label the emergence cohort's deletion topology, and
they disagree. Neither table's own deletion windows can settle it, since those
windows are the quantity under test. The external criterion is an affine-gap
global alignment with free end gaps, recomputed from current UniProt sequences
in r02_align_topology.py.
"""
import csv, json, os
from collections import Counter

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT, U = os.path.join(BASE, "out"), os.environ.get("AFDB_INPUT_DIR", ".")
TOL = 5

rd = lambda p: list(csv.DictReader(open(os.path.join(U, p))))


def classify(start, end, length, tol=TOL):
    n_edge, c_edge = start <= 1 + tol, end >= length - tol
    if n_edge and not c_edge:
        return "N-Terminal Loss"
    if c_edge and not n_edge:
        return "C-Terminal Loss"
    if n_edge and c_edge:
        return "Whole-chain"
    return "Internal Deletion / Splice"


def main():
    s5 = {(r["gene"], r["isoform_acc"]): r for r in rd("emergence_cohort.csv")}
    s7 = {(r["gene"], r["isoform_acc"]): r for r in rd("emergence_deletion_windows.csv")}
    aln = {(r["gene"], r["isoform_acc"]): r
           for r in csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv")))}
    keys = sorted(set(s5) & set(s7) & set(aln))
    if len(keys) != 208:
        raise SystemExit("FATAL: expected 208 joined events, got %d" % len(keys))

    agree5 = agree7 = 0
    per_event = []
    for k in keys:
        r = aln[k]
        ws = [tuple(map(int, w.split("-"))) for w in r["deletion_window"].split(";") if w]
        start, end = min(a for a, _ in ws), max(b for _, b in ws)
        ref = classify(start, end, int(r["canon_len"]))
        agree5 += ref == s5[k]["truncation_topology"]
        agree7 += ref == s7[k]["topology"]
        per_event.append({
            "gene": k[0], "isoform_acc": k[1], "canon_len": r["canon_len"],
            "window_exact_match_rule": "",
            "window_published_alignment": s7[k]["deletion_window"],
            "window_recomputed": r["deletion_window"],
            "n_blocks_recomputed": r["n_deleted_blocks"],
            "junction_substitutions": r["junction_substitutions"],
            "topology_table_s5": s5[k]["truncation_topology"],
            "topology_table_s7": s7[k]["topology"],
            "topology_recomputed": ref,
        })

    verdict = {
        "n_events": len(keys),
        "criterion": "affine-gap global alignment, free end gaps, BLOSUM62, tol=%d" % TOL,
        "agreement_with_recomputed_topology": {
            "table_s5_exact_match_rule": agree5,
            "table_s7_pipeline_alignment": agree7,
        },
        "topology_counts": {
            "table_s5": dict(Counter(s5[k]["truncation_topology"] for k in keys)),
            "table_s7": dict(Counter(s7[k]["topology"] for k in keys)),
            "recomputed": dict(Counter(e["topology_recomputed"] for e in per_event)),
        },
        "note": ("Both assignments under-call C-terminal loss. The exact-match "
                 "windows stop short of the C-terminus wherever splicing leaves novel "
                 "junction residues, which a strict terminal test then reads as internal."),
    }
    json.dump(verdict, open(os.path.join(OUT, "topology_arbitration.json"), "w"), indent=2)
    with open(os.path.join(OUT, "topology_per_event.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_event[0]))
        w.writeheader(); w.writerows(per_event)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
