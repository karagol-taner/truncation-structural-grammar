"""Figure S1. Controls for the feature-deletion rate and for the length-dilution model.

Inputs: auxiliary/random_window_control.json, analysis/pairs_ipsae.csv
Output: figures/FigureS1_controls.{png,pdf,svg}
"""
import json, os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))

RWC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "auxiliary", "random_window_control.json")
PAIRS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analysis", "pairs_ipsae.csv")
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "figures")
OBS, NULL = "#1f77b4", "#bdbdbd"


def main():
    figstyle.apply()
    rwc = json.load(open(RWC))["results"]["alignment windows"]
    df = pd.read_csv(PAIRS)
    df["in_emergence_cohort"] = df["in_emergence_cohort"].astype(str) == "True"

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6), dpi=200)

    # (a) observed against a matched random-window null
    ax = axes[0]
    labels = ["Any curated\nfeature", "Disorder or\ncompositional bias", "Inhibitory or\nregulatory"]
    obs = [rwc["observed_any_feature"], rwc["observed_disorder_or_compbias"],
           rwc["observed_inhibitory"]]
    exp = [rwc["expected_any_feature"], rwc["expected_disorder_or_compbias"],
           rwc["expected_inhibitory"]]
    x = np.arange(len(labels))
    ax.bar(x - 0.2, obs, 0.4, color=OBS, edgecolor="black", linewidth=1.2, label="Observed")
    ax.bar(x + 0.2, exp, 0.4, color=NULL, edgecolor="black", linewidth=1.2,
           label="Random window, same length")
    for xi, (o, e) in enumerate(zip(obs, exp)):
        ax.text(xi - 0.2, o + 1.4, "%.1f%%" % o, ha="center", va="bottom",
                fontsize=9, fontweight="bold")
        ax.text(xi + 0.2, e + 1.4, "%.1f%%" % e, ha="center", va="bottom",
                fontsize=9, fontweight="bold")
        ax.text(xi, max(o, e) + 7.5, "%+.1f pts" % (o - e), ha="center", va="bottom",
                fontsize=9.5, fontweight="bold",
                color="#2ca02c" if (o - e) > 5 else "#7f7f7f")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Deletions overlapping the feature class (%)")
    ax.set_title("Feature deletion against a matched null")
    ax.set_ylim(0, 108)
    ax.legend(loc="upper right", frameon=True, fontsize=8.5)

    # (b) the length-dilution model against a constant
    ax = axes[1]
    em = df[df["in_emergence_cohort"] & (df["iso_len"] > 0)].copy()
    em["pred"] = np.minimum(em["canon_ipTM"] * em["canon_len"] / em["iso_len"], 0.95)
    mae_model = (em["iso_ipTM"] - em["pred"]).abs().mean()
    const = em["iso_ipTM"].mean()
    mae_const = (em["iso_ipTM"] - const).abs().mean()
    ax.scatter(em["pred"], em["iso_ipTM"], s=34, c=OBS, alpha=0.75,
               edgecolor="white", linewidth=0.3)
    lim = [0.4, 1.0]
    ax.plot(lim, lim, color="#525252", ls=":", lw=1.2, label="Perfect prediction")
    ax.axhline(const, color=figstyle.PAL["Ablation"], ls="--", lw=1.6,
               label="Cohort mean, %.3f" % const)
    ax.set_xlim(lim); ax.set_ylim(0.55, 1.0)
    ax.set_xlabel("Predicted isoform ipTM, length-dilution model")
    ax.set_ylabel("Observed isoform ipTM")
    ax.set_title("The length-dilution model on the emergence cohort")
    txt = ("mean absolute error\nlength model   %.4f\ncohort mean    %.4f"
           % (mae_model, mae_const))
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, fontweight="bold", family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                      edgecolor="black", linewidth=1.0, alpha=0.95))
    ax.legend(loc="lower right", frameon=True, fontsize=8.5)

    fig.tight_layout()
    figstyle.save(fig, "FigureS1_controls", OUTDIR)
    print("wrote FigureS1_controls to %s" % OUTDIR)
    print("  length model MAE %.4f, constant-model MAE %.4f, ratio %.2f"
          % (mae_model, mae_const, mae_model / mae_const))


if __name__ == "__main__":
    main()
