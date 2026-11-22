"""Re-derive deletion windows and topology for every pair by alignment.

Replaces the exact prefix/suffix test, which misfiles a C-terminal truncation
as an internal deletion whenever splicing leaves novel junction residues.
Also records overlap with annotated cleavage sites.
"""
import csv, json, os, re, sys, time
from collections import Counter

from Bio.Align import PairwiseAligner, substitution_matrices

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

U = os.environ.get("AFDB_INPUT_DIR", ".")
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "out")
TOL = 5
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 0


def log(m):
    print("%s  %s" % (time.strftime("%H:%M:%S"), m), flush=True)


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


def parse_sites(cell):
    """UniProt feature cell -> list of (start, end, kind)."""
    out = []
    if not cell:
        return out
    for kind, body in re.findall(r"(TRANSIT|SIGNAL|PROPEP)\s+([^;]*)", cell):
        m = re.search(r"(\d+)\.\.(\d+)", body)
        if m:
            out.append((int(m.group(1)), int(m.group(2)), kind))
        elif re.match(r"^\s*\d+\s*$", body):
            p = int(body.strip())
            out.append((p, p, kind))
    return out


def make_aligner():
    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.mode = "global"
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.end_gap_score = 0.0
    return a


def deleted_intervals(aln, canon_len):
    """Canonical (1-based) intervals aligned to gaps in the isoform."""
    tgt, _ = aln.aligned
    covered = [(int(s) + 1, int(e)) for s, e in tgt]
    gaps, prev_end = [], 0
    for s, e in covered:
        if s - 1 > prev_end:
            gaps.append((prev_end + 1, s - 1))
        prev_end = max(prev_end, e)
    if prev_end < canon_len:
        gaps.append((prev_end + 1, canon_len))
    return gaps, covered


def classify(intervals, length, tol):
    if not intervals:
        return "No deletion"
    start = min(s for s, _ in intervals)
    end = max(e for _, e in intervals)
    n_edge, c_edge = start <= 1 + tol, end >= length - tol
    if n_edge and not c_edge:
        return "N-Terminal Loss"
    if c_edge and not n_edge:
        return "C-Terminal Loss"
    if n_edge and c_edge:
        return "Whole-chain"
    return "Internal Deletion / Splice"


def main():
    seqs = load_fasta(os.path.join(DATA, "uniprot_sequences.fasta"))
    sites = {}
    with open(os.path.join(DATA, "uniprot_features.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s = []
            for col in ("Transit peptide", "Signal peptide", "Propeptide"):
                s += parse_sites(r.get(col) or "")
            sites[r["Entry"]] = s
    pairs = list(csv.DictReader(open(os.path.join(U, "pairs_master.csv"))))
    if LIMIT:
        pairs = pairs[:LIMIT]
    log("pairs %d, sequences %d, entries with cleavage sites %d"
        % (len(pairs), len(seqs), sum(1 for v in sites.values() if v)))

    aligner = make_aligner()
    rows, skipped = [], 0
    t0 = time.time()
    for i, p in enumerate(pairs):
        ca, ia = p["canonical_acc"], p["isoform_acc"]
        cs, isq = seqs.get(ca), seqs.get(ia)
        if not cs or not isq:
            skipped += 1
            continue
        aln = aligner.align(cs, isq)[0]
        gaps, covered = deleted_intervals(aln, len(cs))
        total_del = sum(e - s + 1 for s, e in gaps)
        subs = sum(1 for a, b in zip(str(aln[0]), str(aln[1]))
                   if a != "-" and b != "-" and a != b)
        topo = classify(gaps, len(cs), TOL)
        ov = []
        for s, e, kind in sites.get(ca, []):
            if any(not (e < gs or s > ge) for gs, ge in gaps):
                ov.append("%s:%d-%d" % (kind, s, e))
        rows.append({
            "gene": p["gene"], "canonical_acc": ca, "isoform_acc": ia,
            "canon_len": len(cs), "iso_len": len(isq),
            "canon_ipTM": p["canon_ipTM"], "iso_ipTM": p["iso_ipTM"],
            "structural_outcome": p["structural_outcome"],
            "pipeline_len_canon": p["canon_len"], "pipeline_len_iso": p["iso_len"],
            "n_deleted_blocks": len(gaps),
            "deletion_window": ";".join("%d-%d" % g for g in gaps),
            "deletion_length": total_del,
            "junction_substitutions": subs,
            "is_exact_subsequence": isq in cs,
            "topology_alignment": topo,
            "cleavage_site_overlap": ";".join(ov),
            "deletion_is_cleaved_peptide": bool(ov) and len(gaps) == 1 and total_del <= 100,
        })
        if (i + 1) % 500 == 0:
            log("%d/%d  %.1fs" % (i + 1, len(pairs), time.time() - t0))

    if not rows:
        raise SystemExit("FATAL: no pairs aligned")
    fn = os.path.join(OUT, "alignment_topology.csv")
    with open(fn, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    summary = {
        "pairs_aligned": len(rows),
        "pairs_skipped_missing_sequence": skipped,
        "topology_by_alignment": dict(Counter(r["topology_alignment"] for r in rows)),
        "exact_subsequences": sum(r["is_exact_subsequence"] for r in rows),
        "pairs_with_junction_substitutions": sum(r["junction_substitutions"] > 0 for r in rows),
        "multi_block_deletions": sum(r["n_deleted_blocks"] > 1 for r in rows),
        "pairs_overlapping_a_cleavage_site": sum(bool(r["cleavage_site_overlap"]) for r in rows),
        "deletion_is_cleaved_peptide": sum(r["deletion_is_cleaved_peptide"] for r in rows),
        "tolerance": TOL,
    }
    with open(os.path.join(OUT, "alignment_topology_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    log("DONE -> %s" % fn)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
