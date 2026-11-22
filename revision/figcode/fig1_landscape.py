"""Figure 1. The canonical-versus-isoform landscape on both confidence metrics.

Input:  analysis/pairs_ipsae.csv
Output: figures/Figure1_landscape_two_metrics.{png,pdf,svg}
"""
import os, sys
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))

IN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "analysis", "pairs_ipsae.csv")
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figures")
LO, HI = 0.6, 0.8


def classify(c, i):
    if pd.isna(c) or pd.isna(i):
        return "Unscored"
    if c >= HI:
        return "Conservation" if i >= HI else ("Destabilization" if i >= LO else "Ablation")
    if c >= LO:
        return "Mid-tier canonical"
    return "Emergence" if i >= LO else "Baseline"


def panel(ax, df, cn, io, title):
    col = {"Conservation": figstyle.PAL["Conservation"],
           "Destabilization": figstyle.PAL["Destabilization"],
           "Ablation": figstyle.PAL["Ablation"],
           "Emergence": figstyle.PAL["Emergence"],
           "Mid-tier canonical": "#bdbdbd",
           "Baseline": "#e0e0e0"}
    df = df.copy()
    df["cls"] = [classify(c, i) for c, i in zip(df[cn], df[io])]
    for cls in ("Baseline", "Mid-tier canonical", "Conservation", "Destabilization",
                "Ablation", "Emergence"):
        sub = df[df["cls"] == cls]
        if sub.empty:
            continue
        ax.scatter(sub[cn], sub[io], s=14 if cls in ("Baseline", "Mid-tier canonical") else 30,
                   c=col[cls], alpha=0.45 if cls in ("Baseline", "Mid-tier canonical") else 0.85,
                   edgecolor="white", linewidth=0.25,
                   label="%s (n=%d)" % (cls, len(sub)), zorder=2 if cls == "Baseline" else 3)
    for v in (LO, HI):
        ax.axhline(v, color="black", ls="--", lw=1.0, alpha=0.45, zorder=1)
        ax.axvline(v, color="black", ls="--", lw=1.0, alpha=0.45, zorder=1)
    ax.plot([0, 1], [0, 1], color="#525252", ls=":", lw=1.0, alpha=0.7, zorder=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("Canonical homodimer %s" % title)
    ax.set_ylabel("Truncated-isoform homodimer %s" % title)
    ax.set_title(title, color=figstyle.METRIC_COLOR.get(title, "black"))
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, markerscale=1.2)


def main():
    figstyle.apply()
    df = pd.read_csv(IN)
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.6), dpi=200)
    panel(axes[0], df, "canon_ipTM", "iso_ipTM", "ipTM")
    panel(axes[1], df, "canon_ipSAE_min", "iso_ipSAE_min", "ipSAE_min")
    fig.tight_layout()
    figstyle.save(fig, "Figure1_landscape_two_metrics", OUTDIR)
    print("wrote Figure1_landscape_two_metrics to %s" % OUTDIR)


if __name__ == "__main__":
    main()
