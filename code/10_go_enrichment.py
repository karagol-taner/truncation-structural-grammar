"""
10_go_enrichment.py
===================
GO / Reactome / KEGG enrichment of the Interface-Unmasking gene cohort
against a custom proteome-matched background, using g:Profiler.


Run order: 10 of 14
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


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
QUERY_FILE = "unmasking_genes_for_GO.txt"
BG_FULL    = "GO_background_genes.txt"
BG_LOW     = "GO_background_genes_lowtier.txt"
AFDB       = "AFDB_Truncation_Stratified_Analysis_FIXED.csv"

SOURCES          = ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"]
ORG              = "hsapiens"
SIG_THRESH       = 0.05
MIN_TERM_SIZE    = 5
MAX_TERM_SIZE    = 1000
MIN_INTERSECTION = 5

OUT_FULL    = "GO_Enrichment_v2_full.csv"
OUT_LOW     = "GO_Enrichment_v2_lowtier.csv"
OUT_ROBUST  = "GO_Enrichment_v2_robust.csv"
OUT_SUMMARY = "GO_Enrichment_v2_summary.txt"


# ---------------------------------------------------------------------------
def read_genes(path):
    return [g.strip() for g in Path(path).read_text().splitlines() if g.strip()]


def run_enrichment(query, background, label):
    gp = GProfiler(return_dataframe=True, user_agent="KaragolEtAl/1.0")
    res = gp.profile(
        organism=ORG,
        query=query,
        background=background,
        sources=SOURCES,
        significance_threshold_method="g_SCS",
        user_threshold=SIG_THRESH,
        all_results=False,
        ordered=False,
        no_evidences=False,
    )
    if res is None or res.empty:
        return pd.DataFrame()

    res = res[
        (res["term_size"]         >= MIN_TERM_SIZE) &
        (res["term_size"]         <= MAX_TERM_SIZE) &
        (res["intersection_size"] >= MIN_INTERSECTION)
    ].copy()

    from statsmodels.stats.multitest import multipletests
    if not res.empty and "p_value" in res.columns:
        bh = multipletests(res["p_value"], method="fdr_bh")[1]
        res["p_value_bh"] = bh

    res.insert(0, "query_label", label)
    res = res.sort_values("p_value")
    return res


def main():
    # ----- Load inputs
    query_full = read_genes(QUERY_FILE)
    bg_full    = read_genes(BG_FULL)
    bg_low     = read_genes(BG_LOW)

    afdb       = pd.read_csv(AFDB)
    low_genes  = set(afdb.loc[afdb["canon_ipTM"] < 0.5, "gene"].dropna().unique())
    query_low  = sorted(set(query_full) & low_genes)

    print(f"Primary query: {len(query_full)} unmasking genes")
    print(f"Primary background: {len(bg_full)} testable canonicals")
    print(f"Robustness query (Low-tier canonicals only): {len(query_low)} genes")
    print(f"Robustness background (Low-tier canonicals): {len(bg_low)} genes")
    print()

    # ----- Primary
    print("Running PRIMARY enrichment (188 vs 4114) ...")
    full_res = run_enrichment(query_full, bg_full, "primary_188vs4114")
    full_res.to_csv(OUT_FULL, index=False)
    print(f"  {len(full_res)} significant terms after filters -> {OUT_FULL}")

    # ----- Robustness
    print("Running ROBUSTNESS enrichment (Low-tier subquery vs Low-tier bg) ...")
    low_res = run_enrichment(query_low, bg_low, "lowtier_subquery")
    low_res.to_csv(OUT_LOW, index=False)
    print(f"  {len(low_res)} significant terms after filters -> {OUT_LOW}")

    # ----- Robust intersection
    if not full_res.empty and not low_res.empty:
        robust = full_res.merge(
            low_res[["native", "p_value", "p_value_bh", "intersection_size"]],
            on="native", suffixes=("_full", "_low"),
        )
    else:
        robust = pd.DataFrame()
    robust.to_csv(OUT_ROBUST, index=False)
    print(f"  {len(robust)} terms significant in both runs -> {OUT_ROBUST}")

    # ----- Methods-ready summary
    lines = []
    lines.append("ISSUE 2 -- GO ENRICHMENT v2 (g:Profiler) -- SUMMARY")
    lines.append("=" * 60)
    lines.append("")
    lines.append("METHODS BLOCK (paste into Methods)")
    lines.append("-" * 40)
    lines.append(
        "Functional enrichment was performed with g:Profiler "
        "(gprofiler-official Python client) against a custom background "
        "comprising the unique gene symbols of all canonical proteins that "
        f"entered the AFDB-multimer pipeline (n={len(bg_full)}). "
        "Gene Ontology (Biological Process, Molecular Function, Cellular "
        "Component), Reactome and KEGG annotations were queried "
        "simultaneously. Multiple-testing correction used g:SCS at "
        f"alpha={SIG_THRESH}; Benjamini-Hochberg FDR is reported as a "
        "second column. Terms were filtered to annotation size 5 <= |T| <= 1000 "
        "and intersection size >= 5 to exclude uninformative super-terms and "
        "noise. A robustness re-analysis was performed on the strict "
        "gain-of-function subset, restricting both query and background to "
        "canonicals with ipTM < 0.5; only terms significant in BOTH analyses "
        "are reported as 'robust'.")
    lines.append("")
    lines.append("HEADLINE NUMBERS")
    lines.append("-" * 40)
    lines.append(f"Primary enrichment: {len(full_res)} significant terms.")
    lines.append(f"Low-tier robustness: {len(low_res)} significant terms.")
    lines.append(f"Robust (both): {len(robust)} terms.")
    lines.append("")
    if not full_res.empty:
        lines.append("TOP 20 PRIMARY TERMS (by p_value)")
        lines.append("-" * 40)
        cols = ["source", "native", "name", "p_value", "p_value_bh",
                "term_size", "intersection_size"]
        cols = [c for c in cols if c in full_res.columns]
        lines.append(full_res[cols].head(20).to_string(index=False))
        lines.append("")
    if not robust.empty:
        lines.append("TOP 15 ROBUST TERMS (significant in BOTH runs)")
        lines.append("-" * 40)
        cols = ["source", "native", "name", "p_value_full", "p_value_low",
                "intersection_size_full", "intersection_size_low"]
        cols = [c for c in cols if c in robust.columns]
        lines.append(robust[cols].head(15).to_string(index=False))
        lines.append("")
    lines.append("OUTPUTS")
    lines.append("-" * 40)
    lines.append(f"  {OUT_FULL}")
    lines.append(f"  {OUT_LOW}")
    lines.append(f"  {OUT_ROBUST}")

    text = "\n".join(lines)
    Path(OUT_SUMMARY).write_text(text)
    print()
    print(text)


if __name__ == "__main__":
    main()
