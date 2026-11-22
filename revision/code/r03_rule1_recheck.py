"""Test rule 1 (topological asymmetry) on the alignment-derived topology."""
import csv, json, os, re
from collections import Counter, defaultdict

from Bio.Align import PairwiseAligner, substitution_matrices

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT = os.path.join(BASE, "out")
U = os.environ.get("AFDB_INPUT_DIR", ".")
rows = list(csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv"))))
for r in rows:
    r["canon_ipTM"] = float(r["canon_ipTM"]); r["iso_ipTM"] = float(r["iso_ipTM"])
    r["canon_len"] = int(r["canon_len"]); r["deletion_length"] = int(r["deletion_length"])

def win_list(s):
    return [tuple(map(int, w.split("-"))) for w in s.split(";") if w]

def classify(gaps, L, tol):
    if not gaps: return "No deletion"
    s, e = min(a for a, _ in gaps), max(b for _, b in gaps)
    n, c = s <= 1 + tol, e >= L - tol
    if n and not c: return "N-Terminal Loss"
    if c and not n: return "C-Terminal Loss"
    if n and c: return "Whole-chain"
    return "Internal Deletion / Splice"

print("=== proteome-wide topology, tolerance swept (n=%d) ===" % len(rows))
print("tol     N-term    C-term   internal  whole   C-term %")
sweep = []
for tol in (0, 1, 2, 3, 5, 8, 10, 15, 20):
    c = Counter(classify(win_list(r["deletion_window"]), r["canon_len"], tol) for r in rows)
    pct = 100.0 * c["C-Terminal Loss"] / len(rows)
    sweep.append({"tol": tol, **{k: v for k, v in c.items()}, "c_term_pct": round(pct, 1)})
    print("%3d  %8d %9d %10d %6d   %5.1f%%" % (tol, c["N-Terminal Loss"],
          c["C-Terminal Loss"], c["Internal Deletion / Splice"], c["Whole-chain"], pct))
print("\npipeline claim: C-terminal losses are rare, about 6%% of pairs")

# --- rule 1 contingency on high-tier canonicals, as the paper framed it ---
HI = 0.8
hi = [r for r in rows if r["canon_ipTM"] >= HI]
def outcome(r):
    return r["structural_outcome"]
print("\n=== high-tier canonicals (canon_ipTM >= %.1f): n=%d ===" % (HI, len(hi)))
print("pipeline outcome classes:", dict(Counter(outcome(r) for r in hi)))

for tol in (5,):
    tab = defaultdict(Counter)
    for r in hi:
        tab[classify(win_list(r["deletion_window"]), r["canon_len"], tol)][outcome(r)] += 1
    allc = Counter(outcome(r) for r in hi)
    tot = sum(allc.values())
    print("\ntopology x outcome, alignment topology, tol=%d" % tol)
    cons_key = [k for k in allc if "onserv" in k]
    cons = cons_key[0] if cons_key else None
    for topo in ("N-Terminal Loss", "C-Terminal Loss", "Internal Deletion / Splice", "Whole-chain"):
        sub = tab[topo]; n = sum(sub.values())
        if not n: continue
        line = "  %-28s n=%5d" % (topo, n)
        if cons:
            obs = sub[cons] / n
            exp = allc[cons] / tot
            line += "   conservation %5.1f%%  enrichment %.2f-fold" % (100 * obs, obs / exp)
        print(line)
    if cons:
        print("  pipeline enrichment for N-terminal conservation: 1.58-fold")

# --- cleaved targeting peptides inside the emergence cohort ---
s5 = {(r["gene"], r["isoform_acc"]) for r in csv.DictReader(open(os.path.join(U, "emergence_cohort.csv")))}
em = [r for r in rows if (r["gene"], r["isoform_acc"]) in s5]
print("\n=== emergence cohort (%d of 208 matched) ===" % len(em))
cl = [r for r in em if r["cleavage_site_overlap"]]
print("events whose deletion overlaps an annotated cleavage site: %d" % len(cl))
print("  of those, single-block and <=100 residues (the PRDX5/FH shape): %d"
      % sum(1 for r in em if r["deletion_is_cleaved_peptide"] == "True"))
print("  kinds:", dict(Counter(k.split(":")[0] for r in cl for k in r["cleavage_site_overlap"].split(";"))))
print("\ntopology of the emergence cohort by alignment:",
      dict(Counter(classify(win_list(r["deletion_window"]), r["canon_len"], 5) for r in em)))

# --- shared splice events among paralogues in the emergence cohort ---
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.mode = "global"; al.open_gap_score = -11; al.extend_gap_score = -1; al.end_gap_score = 0.0
seqs, acc, buf = {}, None, []
for line in open(os.path.join(BASE, "data", "uniprot_sequences.fasta")):
    if line.startswith(">"):
        if acc: seqs[acc] = "".join(buf)
        acc = line.split("|")[1] if "|" in line else line[1:].split()[0]; buf = []
    elif acc: buf.append(line.strip())
if acc: seqs[acc] = "".join(buf)

groups = defaultdict(list)
for r in em:
    groups[r["deletion_length"]].append(r)
clusters, singles = [], 0
for dl, grp in groups.items():
    if len(grp) == 1:
        singles += 1; continue
    used = set()
    for i, a in enumerate(grp):
        if i in used: continue
        cl_i = [a]; used.add(i)
        for j in range(i + 1, len(grp)):
            if j in used: continue
            b = grp[j]
            sa, sb = seqs.get(a["canonical_acc"]), seqs.get(b["canonical_acc"])
            if not sa or not sb: continue
            aln = al.align(sa, sb)[0]
            ident = sum(1 for x, y in zip(str(aln[0]), str(aln[1])) if x == y and x != "-")
            if ident / max(len(sa), len(sb)) >= 0.90:
                cl_i.append(b); used.add(j)
        if len(cl_i) > 1:
            clusters.append(cl_i)
print("\n=== shared splice events among paralogues (>=90%% canonical identity, same deletion length) ===")
print("clusters found: %d, collapsing %d events into %d" %
      (len(clusters), sum(len(c) for c in clusters), len(clusters)))
for c in sorted(clusters, key=len, reverse=True)[:6]:
    print("   %s  (deletion %d aa)" % (", ".join(x["gene"] for x in c), c[0]["deletion_length"]))
indep = len(em) - sum(len(c) for c in clusters) + len(clusters)
print("independent emergence events after collapsing: %d (pipeline n = 208)" % indep)

json.dump({"topology_sweep": sweep,
           "emergence_cleavage_overlap": len(cl),
           "emergence_independent_events": indep,
           "paralogue_clusters": [[x["gene"] for x in c] for c in clusters]},
          open(os.path.join(OUT, "rule1_recheck.json"), "w"), indent=2)
