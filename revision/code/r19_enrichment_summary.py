"""Summarise the g:Profiler enrichment exports for the two cohorts.

Inputs:  auxiliary/gprofiler/query_*.txt, auxiliary/gprofiler/result_*.csv
Output:  auxiliary/enrichment_summary.json
"""
import csv, json, os

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
GDIR = os.path.join(ROOT, "auxiliary", "gprofiler")

CATALYTIC = "GO:0003824"
HOMODIMER = "GO:0042803"


def genes(name):
    p = os.path.join(GDIR, name)
    return {l.strip() for l in open(p, encoding="utf-8-sig") if l.strip()}


def terms(name):
    p = os.path.join(GDIR, name)
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


def main():
    em = terms("result_emergence_188.csv")
    ab = terms("result_complete_ablation.csv")
    ab_genes = genes("query_complete_ablation.txt")
    rep = genes("query_emergence_reproducing_on_ipSAE.txt")

    def hits(rows, term_id):
        for r in rows:
            if r["term_id"] == term_id:
                return int(r["intersection_size"]), float(r["adjusted_p_value"])
        return None, None

    cat_n, _ = hits(ab, CATALYTIC)
    hom_n, hom_p = hits(ab, HOMODIMER)
    ugt_driven = sum(1 for r in ab
                     if all(g.startswith("UGT") for g in r["intersections"].split(",")))

    out = {
        "enrich_emergence_terms": len(em),
        "enrich_emergence_top_p": round(min(float(r["adjusted_p_value"]) for r in em), 4) if em else None,
        "enrich_emergence_query_size": max(int(r["query_size"]) for r in em) if em else None,
        "enrich_ablation_terms": len(ab),
        "enrich_ablation_min_p": min(float(r["adjusted_p_value"]) for r in ab) if ab else None,
        "enrich_ablation_query_size": max(int(r["query_size"]) for r in ab) if ab else None,
        "enrich_ablation_catalytic_hits": cat_n,
        "enrich_ablation_homodimer_p": round(hom_p, 4) if hom_p is not None else None,
        "enrich_ablation_homodimer_hits": hom_n,
        "enrich_ugt_driven_terms": ugt_driven,
        "enrich_ugt_driver_genes": sum(1 for g in ab_genes if g.startswith("UGT")),
        "enrich_reproducing36_genes": len(rep),
        "enrich_reproducing36_terms": len(terms("result_emergence_reproducing_on_ipSAE.csv")),
    }
    json.dump(out, open(os.path.join(ROOT, "auxiliary", "enrichment_summary.json"), "w"), indent=1)
    for k, v in out.items():
        print("  %-34s %s" % (k, v))


main()
