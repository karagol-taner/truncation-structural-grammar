"""Figure 4. Which emergence events reproduce on the release's own metric.

Inputs: analysis/pairs_ipsae.csv, analysis/plddt_disorder.csv,
        auxiliary/disorder_control.json
Output: figures/Figure4_emergence_discrimination.{png,pdf,svg}
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
PLDDT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analysis", "plddt_disorder.csv")
DISJ = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "auxiliary", "disorder_control.json")
OUTDIR = sys.argv[4] if len(sys.argv) > 4 else os.path.join(ROOT, "figures")
SURV, OTHER = "#1f77b4", "#bdbdbd"


def main():
    figstyle.apply()
    rng = np.random.default_rng(0)
    df = pd.read_csv(PAIRS)
    df["in_emergence_cohort"] = df["in_emergence_cohort"].astype(str) == "True"
    pl = pd.read_csv(PLDDT)
    pl = pl[pl["status"] == "ok"][["uniprotAccession", "fractionPlddtVeryLow"]]
    df = df.merge(pl, left_on="canonical_acc", right_on="uniprotAccession", how="left")
    df = df.rename(columns={"fractionPlddtVeryLow": "disorder"})
    dis = json.load(open(DISJ))

    em = df[df["in_emergence_cohort"]].copy()
    em["reproduces"] = (em["canon_ipSAE_min"] < 0.6) & (em["iso_ipSAE_min"] >= 0.6)
    surv, fell = em[em["reproduces"]], em[~em["reproduces"]]

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 10.4), dpi=200)

    # (a) score change on the two metrics
    ax = axes[0, 0]
    ax.axhline(0, color="black", lw=1.0, alpha=0.5)
    ax.axvline(0, color="black", lw=1.0, alpha=0.5)
    ax.scatter(fell["d_ipTM"], fell["d_ipSAE"], s=34, c=OTHER, edgecolor="white",
               linewidth=0.3, alpha=0.85, label="does not reproduce (n=%d)" % len(fell))
    ax.scatter(surv["d_ipTM"], surv["d_ipSAE"], s=46, c=SURV, edgecolor="black",
               linewidth=0.5, alpha=0.95, label="reproduces on ipSAE (n=%d)" % len(surv))
    lim = [-0.1, 1.0]
    ax.plot(lim, lim, color="#525252", ls=":", lw=1.0, alpha=0.8)
    ax.set_xlabel("Change in ipTM on truncation")
    ax.set_ylabel("Change in ipSAE$_{min}$ on truncation")
    ax.set_title("Score change, emergence cohort")
    ax.legend(loc="upper left", frameon=True, fontsize=8.5)

    # (b) interface residue change
    ax = axes[0, 1]
    data = [fell["d_iface_residues"].dropna(), surv["d_iface_residues"].dropna()]
    bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.8))
    for patch, c in zip(bp["boxes"], (OTHER, SURV)):
        patch.set_facecolor(c); patch.set_edgecolor("black"); patch.set_linewidth(1.2)
    allv = pd.concat(data)
    lo, hi = allv.quantile(0.02), allv.quantile(0.98)
    pad = 0.30 * (hi - lo)
    for i, d in enumerate(data, start=1):
        ax.scatter(rng.normal(i, 0.055, len(d)), d, s=10, c="black", alpha=0.35, zorder=3)
    ax.axhline(0, color="black", ls="--", lw=1.2, alpha=0.6)
    ax.set_ylim(lo - pad, hi + pad)
    n_off = int((allv < lo - pad).sum() + (allv > hi + pad).sum())
    for i, d in enumerate(data, start=1):
        ax.text(i, hi + pad * 0.55, "median %+.0f" % np.median(d), ha="center",
                va="center", fontsize=9.5, fontweight="bold")
    if n_off:
        ax.text(0.99, 0.01, "%d point%s outside axis" % (n_off, "" if n_off == 1 else "s"),
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
                style="italic", fontweight="normal")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Does not\nreproduce", "Reproduces\non ipSAE"])
    ax.set_ylabel("Change in interface residues")
    ax.set_title("Predicted interface size")

    # (c) disorder comparison
    ax = axes[1, 0]
    ab = df[df["in_ablation_cohort"].astype(str) == "True"]["disorder"].dropna()
    ctrl = df[(df["canon_ipTM"] < 0.6) & (~df["in_emergence_cohort"])]["disorder"].dropna()
    emd = em["disorder"].dropna()
    groups = [("Ablation cohort\n(canonical $\\geq$ 0.8)", ab, figstyle.PAL["Ablation"]),
              ("Emergence cohort\n(canonical < 0.6)", emd, figstyle.PAL["Emergence"]),
              ("Matched control\n(canonical < 0.6)", ctrl, "#7f7f7f")]
    bp = ax.boxplot([g[1] for g in groups], widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=1.8))
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(g[2]); patch.set_edgecolor("black"); patch.set_linewidth(1.2)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels([g[0] for g in groups], fontsize=8.5)
    ax.set_ylabel("pLDDT-defined disorder fraction")
    ax.set_title("Canonical disorder by comparison group")
    top = ax.get_ylim()[1]
    ax.plot([1, 1, 2, 2], [top*0.86]*4, lw=1.3, color="black")
    ax.text(1.5, top*0.87, "d = %+.2f" % dis["unmatched"]["cohens_d"], ha="center",
            va="bottom", fontsize=9.5, fontweight="bold")
    ax.plot([2, 2, 3, 3], [top*0.96]*4, lw=1.3, color="black")
    ax.text(2.5, top*0.97, "d = %+.2f" % dis["matched"]["cohens_d"], ha="center",
            va="bottom", fontsize=9.5, fontweight="bold")
    ax.set_ylim(top=top*1.1)

    # (d) stratified by canonical ipTM
    ax = axes[1, 1]
    b = pd.DataFrame(dis["binned"])
    x = np.arange(len(b))
    ax.bar(x - 0.19, b["mean_a"], 0.38, color=figstyle.PAL["Emergence"],
           edgecolor="black", linewidth=1.1, label="Emergence")
    ax.bar(x + 0.19, b["mean_b"], 0.38, color="#7f7f7f",
           edgecolor="black", linewidth=1.1, label="Matched control")
    for xi, row in zip(x, b.itertuples()):
        mark = "d = %+.2f" % row.cohens_d + ("*" if row.p < 0.05 else "")
        ax.text(xi, max(row.mean_a, row.mean_b) + 0.012, mark, ha="center",
                va="bottom", fontsize=8.5, fontweight="bold")
        ax.text(xi - 0.19, 0.004, "%d" % row.n_a, ha="center", va="bottom",
                fontsize=7.5, color="white", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace("ipTM ", "") for n in b["name"]], fontsize=8.5)
    ax.set_xlabel("Canonical ipTM stratum")
    ax.set_ylabel("pLDDT-defined disorder fraction")
    ax.set_title("Disorder within matched confidence strata")
    ax.legend(loc="upper left", frameon=True, fontsize=8.5)
    ax.set_ylim(0, max(b[["mean_a", "mean_b"]].max()) * 1.32)

    fig.tight_layout()
    figstyle.save(fig, "Figure4_emergence_discrimination", OUTDIR)
    print("wrote Figure4_emergence_discrimination to %s" % OUTDIR)


if __name__ == "__main__":
    main()
