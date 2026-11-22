"""Figure 6. Functional enrichment of the emergence and complete-ablation cohorts.

Input:  auxiliary/gprofiler/ (query, background and result files)
Output: figures/Figure6_enrichment_contrast.{png,pdf,svg}
"""
import csv, os, sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))

GDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "auxiliary", "gprofiler")
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figures")

COHORTS = [
    ("query_emergence_208.txt", "result_emergence_188.csv",
     "Emergence", figstyle.PAL["Emergence"]),
    ("query_complete_ablation.txt", "result_complete_ablation.csv",
     "Complete ablation", figstyle.PAL["Ablation"]),
    ("query_emergence_reproducing_on_ipSAE.txt", None,
     "Emergence reproducing\non ipSAE", figstyle.PAL["Baseline"]),
]


def genes(path):
    return {l.strip() for l in open(path, encoding="utf-8") if l.strip()}


def terms(path):
    if path is None or not os.path.exists(path):
        return []
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if r["domain_scope"] if "domain_scope" in r else True:
            out.append({
                "source": r["source"], "name": r["term_name"],
                "p": float(r["adjusted_p_value"]),
                "n_i": int(r["intersection_size"]),
                "inter": set(r["intersections"].split(",")),
            })
    return out


def main():
    figstyle.apply()
    data = []
    for qf, rf, label, colour in COHORTS:
        q = genes(os.path.join(GDIR, qf))
        t = terms(os.path.join(GDIR, rf) if rf else None)
        short = {"Emergence": "Emergence",
                 "Complete ablation": "Complete\nablation"}.get(
                     label, "Emergence\nreproducing\non ipSAE")
        data.append({"label": label, "short": short, "colour": colour,
                     "n_genes": len(q), "terms": t})

    ugt = {g for g in genes(os.path.join(GDIR, COHORTS[1][0]))
           if g.startswith("UGT")}

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.2), dpi=200,
                             gridspec_kw={"width_ratios": [1.62, 1.0]})

    # (a) every significant term, both cohorts, against the significance threshold
    ax = axes[0]
    rows = []
    for d in data:
        for t in sorted(d["terms"], key=lambda x: x["p"], reverse=True):
            rows.append((d["label"], d["colour"], t))
    y = np.arange(len(rows))
    for yi, (lab, col, t) in zip(y, rows):
        x = -np.log10(t["p"])
        driven = t["inter"] <= ugt
        ax.plot([0, x], [yi, yi], color=col, lw=1.6, alpha=0.45, zorder=1)
        ax.scatter([x], [yi], s=26 + 7.0 * t["n_i"], color="white" if driven else col,
                   edgecolor=col, linewidth=1.8, zorder=3,
                   hatch="///" if driven else None)
    ax.axvline(-np.log10(0.05), ls="--", lw=1.3, color="black", zorder=0)
    ax.text(-np.log10(0.05), len(rows) - 0.3, " p = 0.05", fontsize=9,
            fontweight="bold", va="top")
    ax.set_yticks(y)
    ax.set_yticklabels(["%s  (%s)" % (t["name"][:46], t["source"]) for _, _, t in rows],
                       fontsize=8.4)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlabel(r"$-\log_{10}$ adjusted $p$")
    ax.set_title("Significant terms, custom background, annotated domain")
    handles = [plt.Line2D([], [], marker="o", ls="", markersize=9,
                          markerfacecolor=d["colour"], markeredgecolor=d["colour"],
                          label="%s (n = %d genes)" % (d["label"].replace("\n", " "),
                                                       d["n_genes"]))
               for d in data if d["terms"]]
    handles.append(plt.Line2D([], [], marker="o", ls="", markersize=9, markerfacecolor="white",
                              markeredgecolor="black", label="intersection entirely UGT"))
    ax.legend(handles=handles, loc="lower right", frameon=True)

    # (b) cohort size against the number of significant terms
    ax = axes[1]
    xs = np.arange(len(data))
    n_terms = [len(d["terms"]) for d in data]
    ax.bar(xs, n_terms, 0.58, color=[d["colour"] for d in data],
           edgecolor="black", linewidth=1.2)
    for xi, n in zip(xs, n_terms):
        ax.text(xi, n + max(n_terms) * 0.025, "%d" % n, ha="center", va="bottom",
                fontsize=13, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels(["%s\n%d genes" % (d["short"], d["n_genes"]) for d in data],
                       fontsize=8.8)
    ax.set_ylim(0, max(n_terms) * 1.16)
    ax.set_ylabel("Significant terms")
    ax.set_title("Cohort size against recovered enrichment")

    fig.tight_layout()
    figstyle.save(fig, "Figure6_enrichment_contrast", OUTDIR)
    print("terms per cohort:", {d["label"].replace("\n", " "): len(d["terms"]) for d in data})
    print("genes per cohort:", {d["label"].replace("\n", " "): d["n_genes"] for d in data})
    print("UGT genes in the ablation query:", sorted(ugt))


if __name__ == "__main__":
    main()
