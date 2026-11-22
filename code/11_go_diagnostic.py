"""
11_go_diagnostic.py
===================
Diagnostic sensitivity sweep of the unmasking-cohort enrichment.


Run order: 11 of 14
"""

import sys
from pathlib import Path

import pandas as pd

try:
    from gprofiler import GProfiler
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q",
                           "gprofiler-official"])
    from gprofiler import GProfiler


QUERY_FILE = "unmasking_genes_for_GO.txt"
BG_FULL    = "GO_background_genes.txt"

OUT_RELAXED   = "GO_Diagnostic_relaxed.csv"
OUT_BY_SOURCE = "GO_Diagnostic_by_source.csv"
OUT_SUMMARY   = "GO_Diagnostic_summary.txt"

SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"]


def read_genes(p):
    return [g.strip() for g in Path(p).read_text().splitlines() if g.strip()]


def hypergeom_min_hits_for_p(N, K, n, target_p=0.05):
    """
    Minimum number of intersection hits k needed for a hypergeometric
    p-value <= target_p, given background N, term size K, query size n.
    """
    from math import comb
    total = comb(N, n)
    cum = 0.0
    for k in range(min(K, n), -1, -1):
        cum += comb(K, k) * comb(N - K, n - k) / total
        if cum > target_p:
            return k + 1
    return 0


def main():
    query = read_genes(QUERY_FILE)
    bg    = read_genes(BG_FULL)
    print(f"Query: {len(query)} genes; Background: {len(bg)} genes.")

    gp = GProfiler(return_dataframe=True, user_agent="KaragolEtAl-Diag/1.0")

    # ---------- Pass 1: maximally inclusive ----------
    print("\nPass 1: maximally inclusive (all terms, no upper size cap)...")
    res = gp.profile(
        organism="hsapiens",
        query=query,
        background=bg,
        sources=SOURCES,
        significance_threshold_method="g_SCS",
        user_threshold=1.0,
        all_results=True,
        ordered=False,
        no_evidences=False,
    )
    if res is None or res.empty:
        print("  No terms returned.")
        res = pd.DataFrame()
    else:
        res = res.sort_values("p_value")

    res.to_csv(OUT_RELAXED, index=False)
    print(f"  {len(res)} total terms tested -> {OUT_RELAXED}")

    # ---------- Pass 2: per-source top picks ----------
    print("\nPass 2: top 15 per source (p_value < 0.20)...")
    by_src_rows = []
    for src in SOURCES:
        if res.empty:
            continue
        sub = res[(res["source"] == src) & (res["p_value"] < 0.20)].head(15).copy()
        sub["rank_in_source"] = range(1, len(sub) + 1)
        by_src_rows.append(sub)
    by_src = (pd.concat(by_src_rows, ignore_index=True)
              if by_src_rows else pd.DataFrame())
    by_src.to_csv(OUT_BY_SOURCE, index=False)

    # ---------- Power note ----------
    N = len(bg)
    n = len(query)
    power_lines = []
    for K in (20, 50, 100, 200, 500):
        k = hypergeom_min_hits_for_p(N, K, n, 0.05)
        power_lines.append(
            f"  Term size {K}: minimum {k} of {n} query genes needed for "
            f"hypergeometric p<=0.05 (raw, no MTC).")

    # ---------- Summary ----------
    out = []
    out.append("ISSUE 2 -- GO DIAGNOSTIC SENSITIVITY SWEEP -- SUMMARY")
    out.append("=" * 60)
    out.append(f"Query genes: {n}; Background: {N}.")
    out.append("")
    out.append("POWER (raw hypergeometric, before any MTC):")
    out.extend(power_lines)
    out.append("")
    if not res.empty:
        sig05 = (res[res["significant"] == True] if "significant" in res.columns
                 else res[res["p_value"] < 0.05])
        sig10 = res[res["p_value"] < 0.10]
        out.append(f"Terms with raw p < 0.05: {len(sig05)}")
        out.append(f"Terms with raw p < 0.10: {len(sig10)}")
        out.append("")
        out.append("TOP 30 TERMS BY P-VALUE (any source, any size, intersection>=1):")
        cols = ["source", "native", "name", "p_value",
                "term_size", "intersection_size"]
        cols = [c for c in cols if c in res.columns]
        out.append(res[cols].head(30).to_string(index=False))
        out.append("")
        out.append("PER-SOURCE TOP 15 (p < 0.20):")
        out.append(
            by_src[cols + (["rank_in_source"] if "rank_in_source" in by_src.columns else [])]
            .to_string(index=False)
            if not by_src.empty else "  (none)"
        )
    out.append("")
    out.append("INTERPRETATION CHECKLIST")
    out.append("-" * 40)
    out.append("(i) If <= 5 terms have p < 0.05 across all sources at this")
    out.append("    relaxed setting, the unmasking cohort is functionally")
    out.append("    diffuse: the mechanism is structural, not pathway-specific.")
    out.append("    -> Drop the GO panel from Figure 4. Pivot to disorder-")
    out.append("    deletion enrichment + STRING + structural exemplars.")
    out.append("")
    out.append("(ii) If a coherent pathway cluster surfaces (e.g., multiple")
    out.append("     related kinase-signaling, apoptosis, or receptor terms")
    out.append("     all hit p < 0.10 with overlapping gene sets), retain a")
    out.append("     modest GO panel showing those terms with full disclosure")
    out.append("     of the relaxed threshold in the figure caption.")
    out.append("")
    out.append("Outputs:")
    out.append(f"  {OUT_RELAXED}     (every term tested)")
    out.append(f"  {OUT_BY_SOURCE}   (top 15 per source, p<0.20)")

    text = "\n".join(out)
    Path(OUT_SUMMARY).write_text(text)
    print()
    print(text)


if __name__ == "__main__":
    main()
