"""Hypergeometric GO enrichment against the custom background of all canonicals."""
import csv, json, math, os
from collections import Counter, defaultdict

from scipy import stats

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT, DATA = os.path.join(BASE, "out"), os.path.join(BASE, "data")
MIN_TERM, MAX_TERM, MIN_HIT = 5, 1000, 5


def bh(pvals):
    idx = sorted(range(len(pvals)), key=lambda i: pvals[i])
    n = len(pvals)
    out = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(idx), start=1):
        q = min(prev, pvals[i] * n / (n - rank + 1))
        out[i] = prev = q
    return out


def main():
    gene_terms = defaultdict(set)
    term_name = {}
    acc_gene = {}
    with open(os.path.join(DATA, "uniprot_go.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = r["Entry"]
            acc_gene[acc] = (r.get("Gene Names (primary)") or acc).split()[0] if (
                r.get("Gene Names (primary)") or "").strip() else acc
            ids = [t.strip() for t in (r.get("Gene Ontology IDs") or "").split(";") if t.strip()]
            gene_terms[acc] = set(ids)
            for col in ("Gene Ontology (biological process)",
                        "Gene Ontology (molecular function)",
                        "Gene Ontology (cellular component)"):
                for item in (r.get(col) or "").split(";"):
                    item = item.strip()
                    if "[GO:" in item:
                        nm, gid = item.rsplit("[", 1)
                        term_name[gid.strip("] ")] = nm.strip()
    annotated = {a for a, t in gene_terms.items() if t}
    print("canonical accessions with GO annotation: %d of %d" % (len(annotated), len(gene_terms)))

    pairs = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"))))
    background = {r["canonical_acc"] for r in pairs} & annotated
    print("background size: %d" % len(background))

    def cohort(sel):
        return {r["canonical_acc"] for r in pairs if sel(r)} & annotated

    cohorts = {
        "emergence cohort (208 events)":
            cohort(lambda r: r["in_emergence_cohort"] == "True"),
        "emergence reproducing on ipSAE":
            cohort(lambda r: r["in_emergence_cohort"] == "True"
                   and float(r["canon_ipSAE_min"]) < 0.6 <= float(r["iso_ipSAE_min"])),
        "complete ablation cohort":
            cohort(lambda r: r["in_ablation_cohort"] == "True"),
    }

    term_bg = defaultdict(set)
    for a in background:
        for t in gene_terms[a]:
            term_bg[t].add(a)

    results = {}
    for label, query in cohorts.items():
        N, K = len(background), len(query)
        rows = []
        for t, members in term_bg.items():
            n_t = len(members)
            if not (MIN_TERM <= n_t <= MAX_TERM):
                continue
            k = len(members & query)
            if k < MIN_HIT:
                continue
            p = stats.hypergeom.sf(k - 1, N, n_t, K)
            rows.append({"term": t, "name": term_name.get(t, ""), "term_size": n_t,
                         "query_hits": k, "query_size": K, "p": p})
        if rows:
            qs = bh([r["p"] for r in rows])
            for r, q in zip(rows, qs):
                r["p_bh"] = q
            rows.sort(key=lambda r: r["p"])
        sig = [r for r in rows if r.get("p_bh", 1) < 0.05]
        results[label] = {"query_size": K, "terms_tested": len(rows), "significant": len(sig),
                          "top": rows[:8]}
        print("\n=== %s ===" % label)
        print("  query genes with annotation: %d   terms tested: %d   significant at BH 0.05: %d"
              % (K, len(rows), len(sig)))
        for r in rows[:6]:
            print("    %-12s %-46s hits %3d/%4d  p %.2e  p_BH %.3f"
                  % (r["term"], r["name"][:46], r["query_hits"], r["term_size"],
                     r["p"], r.get("p_bh", 1)))

    json.dump(results, open(os.path.join(OUT, "go_enrichment.json"), "w"), indent=2)
    print("\nwritten to out/go_enrichment.json")


if __name__ == "__main__":
    main()
