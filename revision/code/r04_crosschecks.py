"""Outcome-label scheme, alignment windows against the exact-match windows, paralogue inflation."""
import csv, os, re
from collections import Counter, defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT, U = os.path.join(BASE, "out"), os.environ.get("AFDB_INPUT_DIR", ".")
rows = list(csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv"))))
for r in rows:
    r["canon_ipTM"], r["iso_ipTM"] = float(r["canon_ipTM"]), float(r["iso_ipTM"])
    r["canon_len"], r["deletion_length"] = int(r["canon_len"]), int(r["deletion_length"])

print("=== (a) structural_outcome labels across all %d pairs ===" % len(rows))
for k, v in Counter(r["structural_outcome"] for r in rows).most_common():
    sub = [r for r in rows if r["structural_outcome"] == k]
    print("  %-34s %5d   canon_ipTM %.2f-%.2f   iso_ipTM %.2f-%.2f"
          % (k, v, min(r["canon_ipTM"] for r in sub), max(r["canon_ipTM"] for r in sub),
             min(r["iso_ipTM"] for r in sub), max(r["iso_ipTM"] for r in sub)))

print("\n=== (b) alignment windows against the exact-match windows (208 events) ===")
s7 = {(r["gene"], r["isoform_acc"]): r for r in csv.DictReader(open(os.path.join(U, "emergence_deletion_windows.csv")))}
alnw = {(r["gene"], r["isoform_acc"]): r for r in rows}
keys = sorted(set(s7) & set(alnw))
def exact_win(r):
    m = re.findall(r"(\d+)\s*[-–]\s*(\d+)", str(r.get("deletion_window", "")))
    return (int(m[0][0]), int(m[0][1])) if m else None
def aln_win(r):
    ws = [tuple(map(int, w.split("-"))) for w in r["deletion_window"].split(";") if w]
    return (min(a for a, _ in ws), max(b for _, b in ws)) if ws else None
same = diff = 0
lendiff = []
examples = []
for k in keys:
    p, m = exact_win(s7[k]), aln_win(alnw[k])
    if p == m: same += 1
    else:
        diff += 1
        if int(float(s7[k]["canon_len"])) != alnw[k]["canon_len"]:
            lendiff.append(k)
        elif len(examples) < 10:
            examples.append((k, p, m, alnw[k]["canon_len"], alnw[k]["n_deleted_blocks"]))
print("identical windows: %d   different: %d   of which canonical length changed in UniProt: %d"
      % (same, diff, len(lendiff)))
print("  examples where the length is unchanged but the window differs:")
for k, p, m, L, nb in examples:
    print("     %-10s %-12s exact %s   alignment %s   canon_len %d  blocks %d"
          % (k[0], k[1], "%d-%d" % p if p else None, "%d-%d" % m if m else None, L, nb))

print("\n=== (c) paralogue and multi-isoform inflation ===")
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.mode = "global"; al.open_gap_score = -11; al.extend_gap_score = -1; al.end_gap_score = 0.0
seqs, acc, buf = {}, None, []
for line in open(os.path.join(BASE, "data", "uniprot_sequences.fasta")):
    if line.startswith(">"):
        if acc: seqs[acc] = "".join(buf)
        acc = line.split("|")[1] if "|" in line else line[1:].split()[0]; buf = []
    elif acc: buf.append(line.strip())
if acc: seqs[acc] = "".join(buf)

def cluster(subset, label):
    groups = defaultdict(list)
    for r in subset:
        groups[r["deletion_length"]].append(r)
    clusters = []
    for dl, grp in groups.items():
        if len(grp) < 2: continue
        used = set()
        for i, a in enumerate(grp):
            if i in used: continue
            cl, used_i = [a], {i}
            for j in range(i + 1, len(grp)):
                if j in used: continue
                sa, sb = seqs.get(a["canonical_acc"]), seqs.get(grp[j]["canonical_acc"])
                if not sa or not sb or a["canonical_acc"] == grp[j]["canonical_acc"]: continue
                aln = al.align(sa, sb)[0]
                ident = sum(1 for x, y in zip(str(aln[0]), str(aln[1])) if x == y and x != "-")
                if ident / max(len(sa), len(sb)) >= 0.90:
                    cl.append(grp[j]); used_i.add(j)
            used |= used_i
            if len(cl) > 1: clusters.append(cl)
    n_genes = len({r["gene"] for r in subset})
    collapsed = len(subset) - sum(len(c) for c in clusters) + len(clusters)
    print("  %s: %d events, %d distinct genes, %d paralogue clusters, %d events after collapsing paralogues"
          % (label, len(subset), n_genes, len(clusters), collapsed))
    for c in sorted(clusters, key=len, reverse=True)[:4]:
        print("     %s  (deletion %d aa)" % (", ".join(sorted(x["gene"] for x in c)), c[0]["deletion_length"]))
    return clusters

s4 = {(r["gene"], r["isoform_acc"]) for r in csv.DictReader(open(os.path.join(U, "destabilization_leads.csv")))}
s5 = {(r["gene"], r["isoform_acc"]) for r in csv.DictReader(open(os.path.join(U, "emergence_cohort.csv")))}
cluster([r for r in rows if (r["gene"], r["isoform_acc"]) in s4], "destabilization leads (Table S4)")
cluster([r for r in rows if (r["gene"], r["isoform_acc"]) in s5], "emergence cohort  (Table S5)")
