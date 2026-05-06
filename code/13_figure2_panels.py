"""
13_figure2_panels.py
====================
Figure 2 -- residues lost / topology composition / fold enrichment.

Run order: 13 of 14
"""

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import poisson


# ---------------------------------------------------------------------------
# Global style: Bold text & 1:1 aspect ratio
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
    "legend.fontsize":  9,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
})

# ---------------------------------------------------------------------------
# Centralized palette / ordering
# ---------------------------------------------------------------------------
PAL = {
    "Conservation":    "#2ca02c",
    "Destabilization": "#ff7f0e",
    "Ablation":        "#d62728",
    "Emergence":       "#1f77b4",
    "Baseline":        "#7f7f7f",
}

OUTCOME_COLOR = {
    "Conserved High-Confidence Dimer":        "#2ca02c",
    "Conserved Moderate Dimer":               "#90c690",
    "Interface Destabilization (Attenuated)": "#ff7f0e",
    "Gain of Function (Interface Unmasking)": "#1f77b4",
    "Ablation from Moderate":                 "#e07f7f",
    "Complete Interface Ablation":            "#d62728",
}
OUTCOME_SHORT = {
    "Conserved High-Confidence Dimer":        "Conservation (H->H)",
    "Conserved Moderate Dimer":               "Conservation (M->M)",
    "Interface Destabilization (Attenuated)": "Destabilization (H->M)",
    "Gain of Function (Interface Unmasking)": "Emergence",
    "Ablation from Moderate":                 "Ablation (M->L)",
    "Complete Interface Ablation":            "Ablation (H->L)",
}

TOPO_ORDER = ["N-Terminal Loss", "Internal Deletion / Splice", "C-Terminal Loss"]
TOPO_SHORT = {
    "N-Terminal Loss":            "N-term",
    "Internal Deletion / Splice": "Internal",
    "C-Terminal Loss":            "C-term",
}
TOPO_COLOR = {
    "N-Terminal Loss":            "#1f77b4",
    "Internal Deletion / Splice": "#7f7f7f",
    "C-Terminal Loss":            "#d62728",
}


# ---------------------------------------------------------------------------
# Panel 2a -- residues-lost box plot
# ---------------------------------------------------------------------------
def make_figure_2a_residues(table1="Table1_State_Transition_Matrix.csv",
                            fixed="AFDB_Truncation_Stratified_Analysis_FIXED.csv",
                            out="figure_2a_residues", ax=None, save=True):
    fx = pd.read_csv(fixed)
    fx["residues_lost"] = fx["canon_len"] - fx["iso_len"]

    high = fx[fx["canon_ipTM"] >= 0.8].copy()

    def cls(iso):
        if iso >= 0.8: return "Conservation"
        if iso >= 0.5: return "Destabilization"
        return "Ablation"
    high["outcome"] = high["iso_ipTM"].apply(cls)
    high = high[high["residues_lost"] > 0]

    order = ["Conservation", "Destabilization", "Ablation"]

    t1 = pd.read_csv(table1)
    t1_high = t1[t1["canon_tier"] == "High"].set_index("biological_outcome")
    t1_means = {
        "Conservation":    float(t1_high.loc["Conserved High-Confidence Dimer", "avg_residues_lost"]),
        "Destabilization": float(t1_high.loc["Interface Destabilization (Attenuated)", "avg_residues_lost"]),
        "Ablation":        float(t1_high.loc["Complete Interface Ablation", "avg_residues_lost"]),
    }

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    sns.boxplot(
        data=high, x="outcome", y="residues_lost",
        order=order, hue="outcome",
        palette={c: PAL[c] for c in order}, legend=False,
        width=0.55, fliersize=2.5, linewidth=1.2, ax=ax,
    )

    for i, c in enumerate(order):
        ax.scatter(i, t1_means[c], marker="v", s=60, color="black", zorder=5)
        ax.annotate(f"{t1_means[c]:.0f}",
                    xy=(i, t1_means[c]), xytext=(i + 0.30, t1_means[c]),
                    fontsize=9, fontweight="bold", va="center")

    ax.set_xlabel("")
    ax.set_ylabel("Residues lost  (canonical - isoform length)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(bottom=-20)
    ax.set_box_aspect(1)

    if is_standalone and save:
        plt.subplots_adjust(left=0.20, bottom=0.15)
        plt.tight_layout()
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.savefig(f"{out}.svg", bbox_inches="tight")
        plt.close()
        print(f"Wrote {out}.pdf, {out}.svg (high={len(high)})")


# ---------------------------------------------------------------------------
# Panel 2b -- stacked bar of transition composition
# ---------------------------------------------------------------------------
def make_figure_2b_topology(in_csv="Table2_Topology_vs_Biological_Outcome.csv",
                            out="figure_2b_topology", ax=None, save=True):
    df = pd.read_csv(in_csv)
    df_trans = df[df["biological_outcome"] != "Non-Dimer Baseline (Low -> Low)"].copy()

    pivot_table = (df_trans
                   .set_index(["truncation_topology", "biological_outcome"])["count"]
                   .unstack().fillna(0).reindex(TOPO_ORDER))

    pct = pivot_table.div(pivot_table.sum(axis=1), axis=0) * 100

    outcome_order = [
        "Conserved High-Confidence Dimer",
        "Conserved Moderate Dimer",
        "Gain of Function (Interface Unmasking)",
        "Interface Destabilization (Attenuated)",
        "Ablation from Moderate",
        "Complete Interface Ablation",
    ]
    pct = pct[outcome_order]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    x = np.arange(len(TOPO_ORDER))
    bottom = np.zeros(len(TOPO_ORDER))

    for outcome in outcome_order:
        vals = pct[outcome].values
        ax.bar(x, vals, 0.62, bottom=bottom,
               label=OUTCOME_SHORT[outcome],
               color=OUTCOME_COLOR[outcome],
               edgecolor="white", linewidth=1.2)
        bottom += vals

    topo_total  = (df.drop_duplicates("truncation_topology")
                     .set_index("truncation_topology")["topo_total"].to_dict())
    trans_total = pivot_table.sum(axis=1).to_dict()

    xlabels = [f"{TOPO_SHORT[t]} loss\nn = {int(trans_total[t]):,}\nof {topo_total[t]:,}"
               for t in TOPO_ORDER]
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontweight="bold")
    ax.set_ylabel("Percent of transition events  (%)")
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)

    if is_standalone:
        ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5),
                  frameon=False, fontsize=9.5, handlelength=1.2,
                  handletextpad=0.6, prop={'weight': 'bold'})
        if save:
            plt.tight_layout()
            plt.savefig(f"{out}.pdf", bbox_inches="tight")
            plt.savefig(f"{out}.svg", bbox_inches="tight")
            plt.close()
            print(f"Wrote {out}.pdf, {out}.svg")
    else:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0),
                  frameon=False, fontsize=9.5, handlelength=1.2,
                  handletextpad=0.6, prop={'weight': 'bold'})


# ---------------------------------------------------------------------------
# Panel 2c -- forest plot of fold-enrichment
# ---------------------------------------------------------------------------
def make_figure_2c_forest(in_csv="Table2_Topology_vs_Biological_Outcome.csv",
                          out="figure_2c_forest", ax=None, save=True):
    df = pd.read_csv(in_csv)

    grand_total   = df.drop_duplicates("truncation_topology")["topo_total"].sum()
    topo_total    = (df.drop_duplicates("truncation_topology")
                       .set_index("truncation_topology")["topo_total"].to_dict())
    outcome_total = df.groupby("biological_outcome")["count"].sum().to_dict()

    display = [
        ("Conserved High-Confidence Dimer",        "Conservation"),
        ("Interface Destabilization (Attenuated)", "Destabilization"),
        ("Complete Interface Ablation",            "Ablation"),
        ("Gain of Function (Interface Unmasking)", "Emergence"),
    ]

    rows = []
    for outcome_full, outcome_label in display:
        out_tot = outcome_total[outcome_full]
        for topo in TOPO_ORDER:
            sub = df[(df["biological_outcome"] == outcome_full) &
                     (df["truncation_topology"] == topo)]
            obs   = int(sub["count"].iloc[0]) if len(sub) else 0
            t_tot = topo_total[topo]
            expected = t_tot * out_tot / grand_total
            fold     = obs / expected if expected > 0 else 0
            if obs > 0:
                lo, hi = poisson.interval(0.95, obs)
                fold_lo, fold_hi = lo / expected, hi / expected
            else:
                fold_lo = 0.0
                fold_hi = 3.689 / expected
            rows.append(dict(outcome=outcome_label, topo=topo,
                             obs=obs, expected=expected,
                             fold=fold, fold_lo=fold_lo, fold_hi=fold_hi))
    edf = pd.DataFrame(rows)

    gap = 1.0
    y_pos, group_mid = [], []
    cur = 0.0
    for _ in display:
        block_start = cur
        for _ in TOPO_ORDER:
            y_pos.append(cur)
            cur += 1
        group_mid.append((block_start + cur - 1) / 2)
        cur += gap

    if ax is None:
        fig, ax = plt.subplots(figsize=(8.8, 7.8), dpi=200)
        is_standalone = True
    else:
        is_standalone = False

    for i, ((outcome_full, outcome_label), gm) in enumerate(zip(display, group_mid)):
        for j, topo in enumerate(TOPO_ORDER):
            r  = edf.iloc[i * len(TOPO_ORDER) + j]
            yy = y_pos[i * len(TOPO_ORDER) + j]
            ax.errorbar(r["fold"], yy,
                        xerr=[[max(0, r["fold"] - r["fold_lo"])],
                              [max(0, r["fold_hi"] - r["fold"])]],
                        fmt="o", color=TOPO_COLOR[topo],
                        ecolor=TOPO_COLOR[topo],
                        markersize=10.0, capsize=5.0, linewidth=3.5,
                        markeredgecolor="black", markeredgewidth=1.5,
                        zorder=3)
            ax.text(min(r["fold_hi"], 4.6) + 0.20, yy,
                    f"{r['obs']:>4d} / {r['expected']:>4.0f}",
                    fontsize=11.5, fontweight="black", va="center", color="black",
                    family="monospace",
                    path_effects=[pe.withStroke(linewidth=3, foreground="white")])

    for (_, label), gm in zip(display, group_mid):
        ax.text(-0.55, gm, label, ha="right", va="center",
                fontsize=12, fontweight="bold")

    ax.axvline(1, color="black", linestyle="--", linewidth=2.0, alpha=0.8, zorder=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([TOPO_SHORT[t] for _ in display for t in TOPO_ORDER],
                       fontsize=10.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Fold-enrichment  (observed / expected)")
    ax.set_xlim(0, 4.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)

    handles = [mpl.lines.Line2D([0], [0], marker="o", linestyle="",
                                color=TOPO_COLOR[t],
                                markeredgecolor="black", markeredgewidth=1.0,
                                markersize=8.5, label=TOPO_SHORT[t])
               for t in TOPO_ORDER]
    ax.legend(handles=handles, loc="lower right", frameon=False,
              fontsize=10, handletextpad=0.5, title="Topology",
              title_fontsize=11, prop={'weight': 'bold'})

    if is_standalone and save:
        plt.subplots_adjust(left=0.22, right=0.98, top=0.98, bottom=0.10)
        plt.savefig(f"{out}.pdf", bbox_inches="tight")
        plt.savefig(f"{out}.svg", bbox_inches="tight")
        plt.close()
        print(f"Wrote {out}.pdf, {out}.svg")


# ---------------------------------------------------------------------------
# Combined figure
# ---------------------------------------------------------------------------
def make_combined_figure(out="figure_2_combined"):
    """Renders panels a and b on the top row, panel c spanning the bottom."""
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(14, 15), dpi=300)
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            height_ratios=[1, 1.2], wspace=0.35, hspace=0.35)

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    make_figure_2a_residues(ax=axA, save=False)
    make_figure_2b_topology(ax=axB, save=False)
    make_figure_2c_forest(ax=axC,  save=False)

    axA.text(-0.25, 1.05, 'a', transform=axA.transAxes,
             fontsize=24, fontweight='bold', va='top', ha='right')
    axB.text(-0.25, 1.05, 'b', transform=axB.transAxes,
             fontsize=24, fontweight='bold', va='top', ha='right')
    axC.text(-0.25, 1.05, 'c', transform=axC.transAxes,
             fontsize=24, fontweight='bold', va='top', ha='right')

    plt.savefig(f"{out}.png", bbox_inches="tight")
    plt.savefig(f"{out}.svg", bbox_inches="tight")
    plt.close()
    print(f"\nWrote combined figure: {out}.png, {out}.svg")


if __name__ == "__main__":
    make_figure_2a_residues()
    make_figure_2b_topology()
    make_figure_2c_forest()
    make_combined_figure()
    print("\nAll tasks completed successfully.")
