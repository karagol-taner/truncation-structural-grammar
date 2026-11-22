"""Figure 3. Deletion topology, its stability, and its effect on outcome.

Inputs: analysis/alignment_topology.csv, auxiliary/topology_arbitration.json,
        auxiliary/rule1_recheck.json
Output: figures/Figure3_topology.{png,pdf,svg}
"""
import json, os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))

ALN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "analysis", "alignment_topology.csv")
ARB = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "auxiliary", "topology_arbitration.json")
SUM = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "auxiliary", "rule1_recheck.json")
OUTDIR = sys.argv[4] if len(sys.argv) > 4 else os.path.join(ROOT, "figures")
ORDER = ["N-Terminal Loss", "C-Terminal Loss", "Internal Deletion / Splice", "Whole-chain"]


def main():
    figstyle.apply()
    aln = pd.read_csv(ALN)
    arb = json.load(open(ARB))
    rec = json.load(open(SUM))

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0), dpi=200)

    # (a) emergence cohort: three assignments side by side
    ax = axes[0]
    sets = [("Exact-match rule", arb["topology_counts"]["table_s5"]),
            ("Pipeline alignment", arb["topology_counts"]["table_s7"]),
            ("Realigned", arb["topology_counts"]["recomputed"])]
    x = np.arange(len(sets))
    bottom = np.zeros(len(sets))
    for topo in ORDER:
        vals = np.array([s[1].get(topo, 0) for s in sets], dtype=float)
        if not vals.any():
            continue
        ax.bar(x, vals, 0.6, bottom=bottom, color=figstyle.TOPO_COLOR[topo],
               edgecolor="white", linewidth=1.2, label=figstyle.TOPO_SHORT[topo])
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 8:
                ax.text(xi, b + v / 2, "%d" % v, ha="center", va="center",
                        color="white", fontsize=9.5, fontweight="bold")
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels([s[0] for s in sets], fontsize=9)
    ax.set_ylabel("Events in the emergence cohort")
    ax.set_ylim(0, 230)
    ax.set_title("Deletion topology, same %d events" % arb["n_events"])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), frameon=False,
              ncol=4, fontsize=8.5, handlelength=1.3, columnspacing=1.0)

    # (b) proteome-wide topology, swept terminal tolerance
    ax = axes[1]
    sweep = pd.DataFrame(rec["topology_sweep"])
    for topo in ORDER:
        if topo not in sweep.columns:
            continue
        ax.plot(sweep["tol"], sweep[topo], marker="o", ms=4.5, lw=1.8,
                color=figstyle.TOPO_COLOR[topo], label=figstyle.TOPO_SHORT[topo])
    ax.set_xlabel("Terminal tolerance (residues)")
    ax.set_ylabel("Pairs")
    ax.set_title("Proteome-wide topology, tolerance swept")
    ax.legend(loc="center right", frameon=True, fontsize=8.5)
    ax.set_ylim(bottom=0)

    # (c) topology effect conditional on deletion size
    ax = axes[2]
    hi = aln[(aln["canon_ipTM"] >= 0.8)
             & aln["topology_alignment"].isin(["N-Terminal Loss", "C-Terminal Loss"])]
    cut = hi["deletion_length"].median()
    groups, vals, ns, cols = [], [], [], []
    for size_lab, mask in (("< %d aa" % cut, hi["deletion_length"] < cut),
                           ("$\\geq$ %d aa" % cut, hi["deletion_length"] >= cut)):
        for topo in ("N-Terminal Loss", "C-Terminal Loss"):
            sub = hi[mask & (hi["topology_alignment"] == topo)]
            if len(sub) < 3:
                continue
            groups.append("%s\n%s" % (figstyle.TOPO_SHORT[topo], size_lab))
            vals.append(100.0 * (sub["iso_ipTM"] >= 0.8).mean())
            ns.append(len(sub)); cols.append(figstyle.TOPO_COLOR[topo])
    x = np.arange(len(vals))
    bars = ax.bar(x, vals, 0.62, color=cols, edgecolor="black", linewidth=1.2)
    for b, v, nn in zip(bars, vals, ns):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.6, "%.1f%%\nn=%d" % (v, nn),
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8.5)
    ax.set_ylabel("Interface conserved (%)")
    ax.set_title("Topology effect depends on deletion size")
    ax.set_ylim(0, max(vals) * 1.35)

    fig.tight_layout()
    figstyle.save(fig, "Figure3_topology", OUTDIR)
    print("wrote Figure3_topology to %s" % OUTDIR)


if __name__ == "__main__":
    main()
