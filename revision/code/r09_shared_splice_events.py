"""Collapse events that share a splice junction across paralogues.

Six UGT1A entries lose the same 89 residues from canonicals that share only
their C-terminal exons, so whole-chain identity does not group them. Events are
grouped on the identity of the deleted segment itself, with the threshold swept.
"""
import csv, json, os, sys
from collections import Counter, defaultdict

from Bio.Align import PairwiseAligner, substitution_matrices

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

U, OUT, DATA = (os.environ.get("AFDB_INPUT_DIR", "."),
                os.path.join(BASE, "out"), os.path.join(BASE, "data"))
IDENTS = (0.80, 0.85, 0.90, 0.95, 0.99)


def load_fasta(path):
    seqs, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc = line.split("|")[1] if "|" in line else line[1:].split()[0]
            buf = []
        elif acc:
            buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def main():
    seqs = load_fasta(os.path.join(DATA, "uniprot_sequences.fasta"))
    aln = {(r["gene"], r["isoform_acc"]): r
           for r in csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv")))}

    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.mode = "global"
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.end_gap_score = 0.0

    def deleted_segment(acc, window):
        seq = seqs.get(acc, "")
        parts = []
        for w in str(window).split(";"):
            if not w:
                continue
            a, b = (int(v) for v in w.split("-"))
            parts.append(seq[a - 1:b])
        return "".join(parts)

    def identity(sa, sb):
        if not sa or not sb:
            return 0.0
        x = al.align(sa, sb)[0]
        same = sum(1 for p, q in zip(str(x[0]), str(x[1])) if p == q and p != "-")
        return same / max(len(sa), len(sb))

    report = {}
    for table, label in (("ablation_cohort.csv", "complete ablation"),
                         ("destabilization_leads.csv", "partial destabilization"),
                         ("emergence_cohort.csv", "emergence")):
        rows = []
        for r in csv.DictReader(open(os.path.join(U, table))):
            a = aln.get((r["gene"], r["isoform_acc"]))
            if a:
                rows.append({"gene": r["gene"], "canonical_acc": a["canonical_acc"],
                             "isoform_acc": r["isoform_acc"],
                             "deletion_length": int(a["deletion_length"]),
                             "deleted_seq": deleted_segment(a["canonical_acc"],
                                                            a["deletion_window"])})
        groups = defaultdict(list)
        for r in rows:
            groups[r["deletion_length"]].append(r)
        per_threshold = {}
        for thr in IDENTS:
            clusters = []
            for _, grp in groups.items():
                if len(grp) < 2:
                    continue
                used = set()
                for i, a in enumerate(grp):
                    if i in used:
                        continue
                    cl = [a]
                    used.add(i)
                    for j in range(i + 1, len(grp)):
                        if j in used or grp[j]["canonical_acc"] == a["canonical_acc"]:
                            continue
                        if identity(a["deleted_seq"], grp[j]["deleted_seq"]) >= thr:
                            cl.append(grp[j])
                            used.add(j)
                    if len(cl) > 1:
                        clusters.append(cl)
            collapsed = len(rows) - sum(len(c) for c in clusters) + len(clusters)
            per_threshold[thr] = {
                "paralogue_clusters": len(clusters),
                "events_in_clusters": sum(len(c) for c in clusters),
                "independent_events": collapsed,
                "clusters": [{"genes": sorted(x["gene"] for x in c),
                              "deletion_length": c[0]["deletion_length"]} for c in clusters],
            }
        genes = len({r["gene"] for r in rows})
        report[label] = {"events": len(rows), "distinct_genes": genes,
                         "by_threshold": {str(k): v for k, v in per_threshold.items()}}
        print("%s: %d events across %d distinct genes" % (label, len(rows), genes))
        for thr in IDENTS:
            d = per_threshold[thr]
            print("   deleted-segment identity >= %.2f : %2d clusters, %3d independent events"
                  % (thr, d["paralogue_clusters"], d["independent_events"]))
        for c in sorted(per_threshold[0.90]["clusters"], key=lambda c: -len(c["genes"]))[:4]:
            print("      %s  (deletion %d aa)" % (", ".join(c["genes"]), c["deletion_length"]))
    json.dump(report, open(os.path.join(OUT, "shared_splice_events.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
