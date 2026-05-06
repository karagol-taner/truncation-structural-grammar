"""
08_idr_binding_profiling.py
===========================
Comprehensive IDR & binding profiling pipeline (AIUPred + AlphaFold).

Required dependencies (Colab cells):
  !pip install -q metapredict requests pandas numpy matplotlib seaborn scipy
  !pip install pandas numpy matplotlib seaborn scipy requests
  !pip install git+https://github.com/doszilab/AIUPred.git

Run order: 8 of 14
"""

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy import stats

try:
    from aiupred import AIUPred
    print("Initializing AIUPred (loading neural networks)...")
    ai_predictor = AIUPred()
except ImportError:
    print(" [ERROR] aiupred not found. Run:"
          " pip install git+https://github.com/doszilab/AIUPred.git")
    ai_predictor = None


HEADERS = {
    "User-Agent": "KaragolEtAl-IDR-Comprehensive/2.1",
    "Accept":     "text/x-fasta",
}


# ---------------------------------------------------------------------------
# Network helpers
# ---------------------------------------------------------------------------
def fetch_isoform_fasta(iso_acc, retries=3):
    """Fetch the isoform sequence as a single string from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{iso_acc}.fasta"
    for k in range(retries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=25)
            if r.status_code == 200:
                return "".join(line.strip() for line in r.text.splitlines()
                               if not line.startswith(">")) or None
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 ** k)
                continue
            return None
        except requests.RequestException:
            time.sleep(2 ** k)
            return None


def fetch_alphafold_plddt(canonical_acc):
    """
    Fetch the canonical's AFDB monomer PDB and compute the fraction of CA
    residues with pLDDT < 50.  Strips isoform suffixes ("-2") to avoid 404s.
    """
    if not canonical_acc or pd.isna(canonical_acc):
        return float('nan')

    base_acc = str(canonical_acc).split('-')[0].strip()
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{base_acc}"

    try:
        time.sleep(0.15)
        res = requests.get(api_url, timeout=15)
        if res.status_code != 200:
            print(f" [WARN] AlphaFold API returned {res.status_code} for {base_acc}")
            return float('nan')

        data = res.json()
        if not data:
            return float('nan')

        pdb_url = data[0].get('pdbUrl')
        if not pdb_url:
            return float('nan')

        r_pdb = requests.get(pdb_url, timeout=20)
        if r_pdb.status_code != 200:
            return float('nan')

        plddt_scores = []
        for line in r_pdb.text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                plddt_scores.append(float(line[60:66].strip()))

        if not plddt_scores:
            return float('nan')

        plddt_array = np.array(plddt_scores)
        return float(np.sum(plddt_array < 50.0) / len(plddt_array))

    except Exception as e:
        print(f" [WARN] AlphaFold fetch failed for {base_acc}: {e}")
        return float('nan')


# ---------------------------------------------------------------------------
# Scorers
# ---------------------------------------------------------------------------
def score_aiupred_binding(seq):
    """Maximum AIUPred binding (MoRF) probability across the isoform sequence."""
    if not seq or ai_predictor is None:
        return float("nan")
    try:
        binding_scores = ai_predictor.predict_binding(seq)
        if len(binding_scores) == 0:
            return float("nan")
        return float(np.max(binding_scores))
    except Exception as e:
        print(f" [WARN] AIUPred calculation failed: {e}")
        return float("nan")


def compute_statistics(df, metric_col, metric_name):
    """Welch t / Mann-Whitney U / Cohen's d between Gain and Ablation."""
    valid = df.dropna(subset=[metric_col])
    gain = valid.loc[valid["Cohort"] == "Dimer Gain (Unmasking)", metric_col]
    abl  = valid.loc[valid["Cohort"] == "Interface Ablation",     metric_col]

    if len(gain) < 2 or len(abl) < 2:
        return ["Insufficient data for statistical testing."]

    t, p_t = stats.ttest_ind(gain, abl, equal_var=False)
    u, p_u = stats.mannwhitneyu(gain, abl, alternative="two-sided")
    d = ((gain.mean() - abl.mean()) /
         np.sqrt((gain.std() ** 2 + abl.std() ** 2) / 2))

    return [
        f"--- {metric_name} ---",
        f"Welch's t = {t:.3f}, p = {p_t:.3e}",
        f"Mann-Whitney U = {u:.0f}, p = {p_u:.3e}",
        f"Cohen's d = {d:.3f}",
        "",
    ]


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    print("Loading cohorts...")
    try:
        unmasking = pd.read_csv("Dimer_Gain_Cohort.csv")
        ablation  = pd.read_csv("Ablation_Cohort.csv")
    except FileNotFoundError as e:
        print(f" [ERROR] Missing input file: {e}")
        return

    print(f" Unmasking cohort: n = {len(unmasking)}")
    print(f" Ablation cohort: n = {len(ablation)}")

    rows = []
    for cohort_label, dfc in [
        ("Dimer Gain (Unmasking)", unmasking),
        ("Interface Ablation",     ablation),
    ]:
        print(f"\nProcessing {cohort_label} ({len(dfc)} events)...")
        for i, row in dfc.iterrows():
            iso_acc = row["isoform_acc"]
            canonical_acc = row.get("canonical_acc", "")

            seq = fetch_isoform_fasta(iso_acc)
            time.sleep(0.2)

            aiupred_peak = score_aiupred_binding(seq)
            af_plddt_fdr = fetch_alphafold_plddt(canonical_acc)

            rows.append({
                "Cohort":                       cohort_label,
                "Gene":                         row.get("gene", ""),
                "Canonical_Acc":                canonical_acc,
                "Isoform_Acc":                  iso_acc,
                "Sequence_Length":              len(seq) if seq else None,
                "Peak_AIUPred_Binding_Score":   aiupred_peak,
                "AF_Disorder_Fraction_pLDDT50": af_plddt_fdr,
            })

            if (i + 1) % 10 == 0:
                print(f" ...{i+1}/{len(dfc)} processed")

    idr_df = pd.DataFrame(rows)
    idr_df.to_csv("IDR_Comprehensive_Scores.csv", index=False)
    print("\nWrote IDR_Comprehensive_Scores.csv")

    # ---- Statistics
    out_lines = ["=== COMPREHENSIVE IDR & BINDING PROFILING SUMMARY (AIUPred) ===", ""]
    summary = (
        idr_df.groupby("Cohort")[
            ["Peak_AIUPred_Binding_Score", "AF_Disorder_Fraction_pLDDT50"]
        ].agg(["count", "mean", "std", "median"]).round(4)
    )
    out_lines.append(summary.to_string())
    out_lines.append("")
    out_lines.extend(compute_statistics(
        idr_df, "Peak_AIUPred_Binding_Score", "AIUPred Peak Binding Propensity"))
    out_lines.extend(compute_statistics(
        idr_df, "AF_Disorder_Fraction_pLDDT50", "AlphaFold pLDDT < 50 Fraction"))

    summary_text = "\n".join(out_lines)
    print("\n" + summary_text)
    Path("IDR_Profiling_summary.txt").write_text(summary_text)

    # ---- Plotting
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), dpi=200)
    palette = {
        "Dimer Gain (Unmasking)": "#1f77b4",
        "Interface Ablation":     "#d62728",
    }

    # Panel A -- AIUPred binding
    valid_ai = idr_df.dropna(subset=["Peak_AIUPred_Binding_Score"])
    if not valid_ai.empty:
        sns.violinplot(data=valid_ai, x="Cohort", y="Peak_AIUPred_Binding_Score",
                       hue="Cohort", palette=palette, legend=False,
                       inner="box", linewidth=0.7, ax=axes[0], cut=0)
        sns.stripplot(data=valid_ai, x="Cohort", y="Peak_AIUPred_Binding_Score",
                      size=2.5, alpha=0.45, color="black", jitter=True, ax=axes[0])
        axes[0].axhline(0.5, color="gray", linestyle="--", alpha=0.55, linewidth=0.6)
        axes[0].set_xlabel("")
        axes[0].set_ylabel("Peak AIUPred Binding Propensity")
        axes[0].spines["top"].set_visible(False)
        axes[0].spines["right"].set_visible(False)
    else:
        axes[0].text(0.5, 0.5, "AIUPred\ndata not found", ha='center', va='center')

    # Panel B -- AlphaFold pLDDT
    valid_af = idr_df.dropna(subset=["AF_Disorder_Fraction_pLDDT50"])
    if not valid_af.empty:
        sns.violinplot(data=valid_af, x="Cohort", y="AF_Disorder_Fraction_pLDDT50",
                       hue="Cohort", palette=palette, legend=False,
                       inner="box", linewidth=0.7, ax=axes[1], cut=0)
        sns.stripplot(data=valid_af, x="Cohort", y="AF_Disorder_Fraction_pLDDT50",
                      size=2.5, alpha=0.45, color="black", jitter=True, ax=axes[1])
        axes[1].set_xlabel("")
        axes[1].set_ylabel("Fraction of Residues with pLDDT < 50")
        axes[1].spines["top"].set_visible(False)
        axes[1].spines["right"].set_visible(False)
    else:
        axes[1].text(0.5, 0.5, "AlphaFold data\nnot found", ha='center', va='center')

    plt.tight_layout()
    plt.savefig("IDR_Comprehensive_Profiling.pdf", bbox_inches="tight")
    plt.savefig("IDR_Comprehensive_Profiling.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Wrote IDR_Comprehensive_Profiling.pdf, IDR_Comprehensive_Profiling.png")


if __name__ == "__main__":
    main()
