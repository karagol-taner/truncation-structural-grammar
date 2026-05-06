"""
14_figure3_panels.py
====================
Figure 3 -- truncation topology, structural disorder, binding propensity,
event-level summary and feature-class breakdown for the Emergence cohort.

Five panels plus a combined 2x3 layout:
  3a  Topology distribution: cohort vs proteome-wide reference.
  3b  AlphaFold structural flexibility (fraction of pLDDT < 50) per cohort.
  3c  AIUPred peak MoRF binding propensity per cohort.
  3d  Event-level summary: percent of events deleting any curated
      feature / a disorder / compositional-bias feature / an
      inhibitory-regulatory feature.
  3e  Top deleted UniProt feature classes, split by disorder annotation.
  combined: 2x3 grid, panel e spans the bottom right two columns.

Reads:
  Domain_Annotations_v2.csv
  Domain_Annotations_v2_per_event.csv
  Table2_Topology_vs_Biological_Outcome.csv
  IDR_Comprehensive_Scores.csv

Run order: 14 of 14
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# ---------------------------------------------------------------------------
# Global style: bold text & 1:1 aspect ratio
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family":      "sans-serif",
    "font.weight":      "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "font.size":        9.5,
    "axes.linewidth":   1.2,
    "axes.labelsize":   10.5,
    "xtick.labelsize":  9.0,
    "ytick.labelsize":  9.0,
    "legend.fontsize":  9.0,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
})


# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------
GAIN_PURPLE    = "#6a3d9a"
ABLATION_RED   = "#cb181d"

EMERGENCE_BLUE = "#1f77b4"
BASELINE_GRAY  = "#bbbbbb"
DISORDER_GREEN = "#2ca02c"
NEUTRAL_GRAY   = "#bbbbbb"
INHIB_RED      = "#d62728"
ANY_BLUE       = "#1f77b4"

TOPO_ORDER = ["N-Terminal Loss", "Internal Deletion / Splice", "C-Terminal Loss"]
TOPO_SHORT = {
    "N-Terminal Loss":            "N-terminal",
    "Internal Deletion / Splice": "Internal /\nsplice",
    "C-Terminal Loss":            "C-terminal",
}


# ---------------------------------------------------------------------------
# Panel 3a -- Topology distribution
# ---------------------------------------------------------------------------
def make_figure_3a_topology(per_event="Domain_Annotations_v2_per_event.csv",
                            table2="Table2_Topology_vs_Biological_Outcome.csv",
                            out="figure_3a_topology", ax=None, save=True):
    df = pd.read_csv(per_event)
    cohort_counts = df["topology"].value_counts()
    cohort_total  = int(cohort_counts.sum())
    cohort_pct    = (cohort_counts / cohort_total * 100).to_dict()

    t2       = pd.read_csv(table2)
    proteome = (t2.drop_duplicates("truncation_topology")
                  .set_index("truncation_topology")["topo_total"].to_dict())
    proteome_total = int(sum(proteome.values()))
    proteome_pct   = {k: v / proteome_total * 100 for k, v in proteome.items()}

    cohort_vals   = [cohort_pct.get(t, 0)   for t in TOPO_ORDER]
    proteome_vals = [proteome_pct.get(t, 0) for t in TOPO_ORDER]
    cohort_n      = [int(cohort_counts.get(t, 0)) for t in TOPO_ORDER]
    proteome_n    = [int(proteome.get(t, 0))      for t in TOPO_ORDER]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    x     = np.arange(len(TOPO_ORDER))
    width = 0.36

    ax.bar(x - width/2, cohort_vals, width,
           color=EMERGENCE_BLUE, edgecolor="black", linewidth=1.2,
           label=f"Emergence cohort (n={cohort_total})")
    ax.bar(x + width/2, proteome_vals, width,
           color=BASELINE_GRAY, edgecolor="black", linewidth=1.2,
           label=f"Proteome-wide (n={proteome_total:,})")

    for xi, cv, cn in zip(x, cohort_vals, cohort_n):
        ax.text(xi - width/2, cv + 1.5, f"{cv:.1f}%",
                ha="center", va="bottom", fontsize=9.5,
                fontweight="bold", color="#11436b")
        ax.text(xi - width/2, -3.5, f"n={cn}",
                ha="center", va="top", fontsize=8.5,
                fontweight="bold", color="#11436b")

    for xi, pv, pn in zip(x, proteome_vals, proteome_n):
        ax.text(xi + width/2, pv + 1.5, f"{pv:.1f}%",
                ha="center", va="bottom", fontsize=9.5,
                fontweight="bold", color="#555555")

    ax.set_xticks(x)
    ax.set_xticklabels([TOPO_SHORT[t] for t in TOPO_ORDER], fontweight="bold")
    ax.set_ylabel("Percent of events (%)", fontweight="bold")
    ax.set_ylim(0, 75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title('Truncation Topology', fontweight="bold", pad=10)
    ax.set_box_aspect(1)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5,
              prop={'weight': 'bold'})

    if is_standalone and save:
        plt.tight_layout()
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.close()


# ---------------------------------------------------------------------------
# Panel 3b -- AlphaFold pLDDT structural confidence
# ---------------------------------------------------------------------------
def make_figure_3b_alphafold(in_csv="IDR_Comprehensive_Scores.csv",
                             out="figure_3b_alphafold", ax=None, save=True):
    if not os.path.exists(in_csv):
        raise FileNotFoundError(f"\n\nCRITICAL ERROR: '{in_csv}' not found. "
                                "Run profiling pipeline first.")

    idr_data = pd.read_csv(in_csv).dropna(subset=["AF_Disorder_Fraction_pLDDT50"])

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 6), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    sns.violinplot(data=idr_data, x='Cohort', y='AF_Disorder_Fraction_pLDDT50',
                   hue='Cohort',
                   palette=[GAIN_PURPLE, ABLATION_RED], legend=False,
                   inner="box", ax=ax, linewidth=1.2, cut=0)
    sns.stripplot(data=idr_data, x="Cohort", y="AF_Disorder_Fraction_pLDDT50",
                  size=2.5, alpha=0.45, color="black", jitter=True, ax=ax)

    ax.set_title('AlphaFold Structural Flexibility', fontweight="bold", pad=10)
    ax.set_ylabel('Fraction of residues with pLDDT < 50', fontweight="bold")
    ax.set_xlabel('', fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)

    if is_standalone and save:
        plt.tight_layout()
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.close()


# ---------------------------------------------------------------------------
# Panel 3c -- AIUPred peak binding propensity (boxplot)
# ---------------------------------------------------------------------------
def make_figure_3c_aiupred(in_csv="IDR_Comprehensive_Scores.csv",
                           out="figure_3c_aiupred", ax=None, save=True):
    if not os.path.exists(in_csv):
        raise FileNotFoundError(f"\n\nCRITICAL ERROR: '{in_csv}' not found. "
                                "Run profiling pipeline first.")

    idr_data = pd.read_csv(in_csv).dropna(subset=["Peak_AIUPred_Binding_Score"])

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 6), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    sns.boxplot(data=idr_data, x='Cohort', y='Peak_AIUPred_Binding_Score',
                hue='Cohort',
                palette=[GAIN_PURPLE, ABLATION_RED], legend=False,
                showfliers=False, ax=ax, linewidth=1.5)
    sns.stripplot(data=idr_data, x="Cohort", y="Peak_AIUPred_Binding_Score",
                  size=3.0, alpha=0.55, color="black", jitter=True, ax=ax)

    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.6, linewidth=1.2)
    ax.set_title('AIUPred Binding Propensity', fontweight="bold", pad=10)
    ax.set_ylabel('Peak MoRF Binding Probability', fontweight="bold")
    ax.set_xlabel('', fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)

    if is_standalone and save:
        plt.tight_layout()
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.close()


# ---------------------------------------------------------------------------
# Panel 3d -- Event-level summary percentages
# ---------------------------------------------------------------------------
def make_figure_3d_event_summary(per_event="Domain_Annotations_v2_per_event.csv",
                                 out="figure_3d_event_summary", ax=None, save=True):
    pe = pd.read_csv(per_event)

    n_total    = len(pe)
    n_any      = int((pe["n_features_overlapping"] >= 1).sum())
    n_disorder = int(pe["any_disorder_lexical"].sum())
    n_inhib    = int(pe["any_inhibitory_lexical"].sum())

    pct_any      = 100 * n_any      / n_total if n_total else 0
    pct_disorder = 100 * n_disorder / n_total if n_total else 0
    pct_inhib    = 100 * n_inhib    / n_total if n_total else 0

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    cats = [
        "Any curated\nfeature deleted",
        "Disordered /\ncomp.-bias region\ndeleted",
        "Inhibitory /\nregulatory feature\ndeleted",
    ]
    pcts   = [pct_any,  pct_disorder, pct_inhib]
    counts = [n_any,    n_disorder,   n_inhib]
    cols   = [ANY_BLUE, DISORDER_GREEN, INHIB_RED]

    bars = ax.bar(np.arange(3), pcts, color=cols,
                  edgecolor="black", linewidth=1.2, width=0.62)

    for bar, pct, n in zip(bars, pcts, counts):
        ax.text(bar.get_x() + bar.get_width()/2, pct + 2.5, f"{pct:.1f}%",
                ha="center", va="bottom", fontsize=10.5, fontweight="bold")
        ax.text(bar.get_x() + bar.get_width()/2, pct - 5.0, f"{n}/{n_total}",
                ha="center", va="top", fontsize=9.0, fontweight="bold",
                color="white" if pct > 15 else "black")

    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(cats, fontsize=9.0, fontweight="bold")
    ax.set_ylabel("Percent of emergence events (%)", fontweight="bold")
    ax.set_ylim(0, 105)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Event-Level Summary", fontsize=11.0, fontweight="bold", pad=12)
    ax.set_box_aspect(1)

    if is_standalone and save:
        plt.tight_layout()
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.close()


# ---------------------------------------------------------------------------
# Panel 3e -- Feature-class counts
# ---------------------------------------------------------------------------
def make_figure_3e_features(long_csv="Domain_Annotations_v2.csv",
                            out="figure_3e_features", ax=None, save=True):
    lf = pd.read_csv(long_csv)

    split = (lf.groupby(["feature_type", "is_disorder_lexical"])
                .size().unstack(fill_value=0))
    if True  not in split.columns: split[True]  = 0
    if False not in split.columns: split[False] = 0
    split["total"] = split[True] + split[False]
    split = split.sort_values("total", ascending=True)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8.5, 5), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    y    = np.arange(len(split))
    dis  = split[True].values
    nond = split[False].values

    ax.barh(y, dis, color=DISORDER_GREEN, edgecolor="black", linewidth=1.0,
            height=0.70, label="Disorder / compositional-bias\nannotation")
    ax.barh(y, nond, left=dis, color=NEUTRAL_GRAY,
            edgecolor="black", linewidth=1.0, height=0.70,
            label="Other annotation")

    max_val = max(split["total"]) if len(split) > 0 else 10
    ax.set_xlim(0, max_val * 1.35)

    for i, (d, nd) in enumerate(zip(dis, nond)):
        total = int(d + nd)
        ax.text(total + (max_val * 0.02), i, f"{total}",
                va="center", ha="left", fontsize=9.5, fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels(split.index, fontsize=9.5, fontweight="bold")
    ax.set_xlabel(f"Number of intersection events "
                  f"(deletion-window x feature; n={len(lf):,} total)",
                  fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=9.0,
              prop={'weight': 'bold'})
    ax.set_title("Top Deleted UniProt Feature Classes",
                 fontsize=11.0, fontweight="bold", pad=12)

    if is_standalone and save:
        plt.subplots_adjust(left=0.35)
        plt.tight_layout()
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.close()


# ---------------------------------------------------------------------------
# Combined Figure 3 (2x3 grid; panel e spans bottom-right two cols)
# ---------------------------------------------------------------------------
def make_combined_figure_3(out="figure_3_combined"):
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(18, 12), dpi=300)
    gs  = gridspec.GridSpec(2, 3, figure=fig, wspace=0.35, hspace=0.45)

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])
    axD = fig.add_subplot(gs[1, 0])
    axE = fig.add_subplot(gs[1, 1:3])

    make_figure_3a_topology(ax=axA,      save=False)
    make_figure_3b_alphafold(ax=axB,     save=False)
    make_figure_3c_aiupred(ax=axC,       save=False)
    make_figure_3d_event_summary(ax=axD, save=False)
    make_figure_3e_features(ax=axE,      save=False)

    for ax_panel, letter in zip([axA, axB, axC, axD, axE],
                                ['a', 'b', 'c', 'd', 'e']):
        ax_panel.annotate(letter, xy=(-0.15, 1.05), xycoords='axes fraction',
                          fontsize=24, fontweight='bold', va='bottom', ha='right')

    fig.tight_layout(pad=3.0, w_pad=2.5, h_pad=3.0)

    plt.savefig(f"{out}.png", bbox_inches="tight")
    plt.savefig(f"{out}.svg", bbox_inches="tight")
    plt.close()
    print(f"\nWrote combined figure: {out}.png, {out}.svg")


if __name__ == "__main__":
    if not os.path.exists("IDR_Comprehensive_Scores.csv"):
        raise FileNotFoundError(
            "\n\nCRITICAL ERROR: 'IDR_Comprehensive_Scores.csv' not found.\n"
            "Strict policy prohibits the use of placeholder or generated mock data.\n"
        )

    make_figure_3a_topology()
    make_figure_3b_alphafold()
    make_figure_3c_aiupred()
    make_figure_3d_event_summary()
    make_figure_3e_features()
    make_combined_figure_3()

    print("\nAll visualization tasks completed successfully.")
