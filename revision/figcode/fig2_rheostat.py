"""Figure 2. The graded response measured on interface size.

Inputs: analysis/pairs_ipsae.csv, analysis/alignment_topology.csv
Output: figures/Figure2_graded_interface_loss.{png,pdf,svg}
"""
import json, os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))

PAIRS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "analysis", "pairs_ipsae.csv")
ALN = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analysis", "alignment_topology.csv")
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "figures")


def main():
    figstyle.apply()
    df = pd.read_csv(PAIRS)
    aln = pd.read_csv(ALN)[["gene", "isoform_acc", "topology_alignment",
                            "deletion_length", "canon_len"]]
    df = df.merge(aln, on=["gene", "isoform_acc"], how="left", suffixes=("", "_a"))
    df = df[(df["canon_ipSAE_min"] >= 0.6) & (df["canon_iface_residues"] > 0)
            & df["iso_iface_residues"].notna() & (df["canon_len_a"] > 0)].copy()
    df["frac_deleted"] = df["deletion_length"] / df["canon_len_a"]
    df["iface_ratio"] = df["iso_iface_residues"] / df["canon_iface_residues"]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.2), dpi=200)

    # (a) interface retained against how much of the chain is removed
    ax = axes[0]
    bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 1.01)]
    labels, data = [], []
    for lo, hi in bins:
        sub = df[(df["frac_deleted"] >= lo) & (df["frac_deleted"] < hi)]["iface_ratio"]
        if len(sub) < 5:
            continue
        labels.append("%.0f-%.0f%%\nn=%d" % (100 * lo, 100 * min(hi, 1.0), len(sub)))
        data.append(sub.values)
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.8))
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(data)))
    for patch, c in zip(bp["boxes"], cmap):
        patch.set_facecolor(c); patch.set_edgecolor("black"); patch.set_linewidth(1.1)
    meds = [np.median(d) for d in data]
    ax.plot(np.arange(1, len(meds) + 1), meds, color="black", lw=1.6, ls="--", zorder=4)
    for i, m in enumerate(meds, start=1):
        ax.text(i, m + 0.045, "%.3f" % m, ha="center", va="bottom",
                fontsize=9, fontweight="bold")
    ax.set_xticks(np.arange(1, len(labels) + 1)); ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Chain removed by truncation")
    ax.set_ylabel("Interface residues retained (ratio)")
    ax.set_title("Graded loss of predicted interface")
    ax.set_ylim(-0.05, 1.18)

    # (b) confidence change against interface-size change
    ax = axes[1]
    ax.scatter(df["d_iface_residues"], df["d_ipTM"], s=16,
               c=figstyle.METRIC_COLOR["ipTM"], alpha=0.45, edgecolor="none",
               label="ipTM")
    ax.scatter(df["d_iface_residues"], df["d_ipSAE"], s=16,
               c=figstyle.METRIC_COLOR["ipSAE_min"], alpha=0.45, edgecolor="none",
               label="ipSAE$_{min}$")
    ax.axhline(0, color="black", lw=1.0, alpha=0.5)
    ax.axvline(0, color="black", lw=1.0, alpha=0.5)
    ax.set_xlabel("Change in interface residues")
    ax.set_ylabel("Change in confidence score")
    ax.set_title("Confidence tracks interface size")
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    ranks = df[["d_iface_residues", "d_ipTM", "d_ipSAE"]].rank()
    rho_iptm = ranks["d_iface_residues"].corr(ranks["d_ipTM"])
    rho_ipsae = ranks["d_iface_residues"].corr(ranks["d_ipSAE"])
    ax.text(0.03, 0.97, "Spearman\nipTM   %+.3f\nipSAE  %+.3f" % (rho_iptm, rho_ipsae),
            transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
            fontweight="bold", family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="black", lw=1.0, alpha=0.95))

    fig.tight_layout()
    figstyle.save(fig, "Figure2_graded_interface_loss", OUTDIR)
    print("wrote Figure2_graded_interface_loss to %s" % OUTDIR)


if __name__ == "__main__":
    main()
