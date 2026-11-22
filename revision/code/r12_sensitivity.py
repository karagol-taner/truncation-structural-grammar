"""Sensitivity of the headline results to the artifact classes, with bootstrap intervals."""
import csv, json, math, os, random
from collections import Counter

from scipy import stats

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT = os.path.join(BASE, "out")
BOOT = 4000
SEED = 20261119


def f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")


def boot_ci(vals, fn, n=BOOT, seed=SEED):
    rnd = random.Random(seed)
    k = len(vals)
    if k < 3:
        return (float("nan"), float("nan"))
    stats_ = []
    for _ in range(n):
        samp = [vals[rnd.randrange(k)] for _ in range(k)]
        stats_.append(fn(samp))
    stats_.sort()
    return (stats_[int(0.025 * n)], stats_[int(0.975 * n)])


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"))))
    aln = {(r["gene"], r["isoform_acc"]): r
           for r in csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv")))}
    for r in rows:
        for k in ("canon_ipTM", "iso_ipTM", "canon_ipSAE_min", "iso_ipSAE_min"):
            r[k] = f(r[k])
        a = aln.get((r["gene"], r["isoform_acc"]))
        r["topology"] = a["topology_alignment"] if a else ""
        r["cleaved"] = bool(a and a["cleavage_site_overlap"])
        r["cleaved_only"] = bool(a and a["deletion_is_cleaved_peptide"] == "True")
    out = {}

    print("=== sensitivity of the headline numbers to excluding the cleaved-peptide class ===")
    subsets = [("all pairs", rows),
               ("excluding any cleavage-site overlap",
                [r for r in rows if not r["cleaved"]]),
               ("excluding single-block cleaved peptides only",
                [r for r in rows if not r["cleaved_only"]])]
    for label, sub in subsets:
        topo = Counter(r["topology"] for r in sub if r["topology"])
        tot = sum(topo.values())
        cterm = 100.0 * topo["C-Terminal Loss"] / tot if tot else float("nan")
        hi = [r for r in sub if r["canon_ipTM"] >= 0.8]
        cons = sum(1 for r in hi if r["iso_ipTM"] >= 0.8)
        em = [r for r in sub if r["canon_ipTM"] < 0.6 <= r["iso_ipTM"]]
        surv = sum(1 for r in em if r["canon_ipSAE_min"] < 0.6 <= r["iso_ipSAE_min"])
        print("  %-44s n=%4d  C-term %.1f%%  high-tier %3d  conserved %.1f%%  emergence %3d  reproducing %3d"
              % (label, len(sub), cterm, len(hi), 100.0*cons/len(hi) if hi else float("nan"),
                 len(em), surv))
        out[label] = {"n": len(sub), "c_terminal_pct": cterm, "high_tier": len(hi),
                      "conserved_pct": 100.0*cons/len(hi) if hi else None,
                      "emergence": len(em), "emergence_reproducing": surv}

    print("\n=== bootstrap intervals on the headline proportions, all pairs ===")
    topo_flags = [1 if r["topology"] == "C-Terminal Loss" else 0 for r in rows if r["topology"]]
    lo, hi_ = boot_ci(topo_flags, lambda v: 100.0 * sum(v) / len(v))
    print("  C-terminal share of pairs:            %.1f%%  95%% CI %.1f to %.1f"
          % (100.0 * sum(topo_flags) / len(topo_flags), lo, hi_))
    out["ci_c_terminal"] = [lo, hi_]

    hi_rows = [r for r in rows if r["canon_ipTM"] >= 0.8]
    for topo in ("N-Terminal Loss", "C-Terminal Loss"):
        flags = [1 if r["iso_ipTM"] >= 0.8 else 0 for r in hi_rows if r["topology"] == topo]
        if len(flags) < 5: continue
        lo, hi_ = boot_ci(flags, lambda v: 100.0 * sum(v) / len(v))
        print("  conservation after %-16s %.1f%%  95%% CI %.1f to %.1f  (n=%d)"
              % (topo + ":", 100.0 * sum(flags) / len(flags), lo, hi_, len(flags)))
        out["ci_conservation_%s" % topo.split("-")[0]] = [lo, hi_]

    em = [r for r in rows if r["canon_ipTM"] < 0.6 <= r["iso_ipTM"]]
    flags = [1 if r["canon_ipSAE_min"] < 0.6 <= r["iso_ipSAE_min"] else 0 for r in em]
    lo, hi_ = boot_ci(flags, lambda v: 100.0 * sum(v) / len(v))
    print("  emergence reproducing on ipSAE:       %.1f%%  95%% CI %.1f to %.1f  (n=%d)"
          % (100.0 * sum(flags) / len(flags), lo, hi_, len(flags)))
    out["ci_emergence_reproducing"] = [lo, hi_]

    print("\n=== is the topology effect independent of deletion size? ===")
    hi_rows2 = [r for r in hi_rows if r["topology"] in ("N-Terminal Loss", "C-Terminal Loss")]
    sizes = sorted(float(aln[(r["gene"], r["isoform_acc"])]["deletion_length"]) for r in hi_rows2)
    cut = sizes[len(sizes) // 2]
    for lab, pick in (("deletion below the median (%d aa)" % cut,
                       lambda r: float(aln[(r["gene"], r["isoform_acc"])]["deletion_length"]) < cut),
                      ("deletion at or above the median",
                       lambda r: float(aln[(r["gene"], r["isoform_acc"])]["deletion_length"]) >= cut)):
        sub = [r for r in hi_rows2 if pick(r)]
        line = "  %-36s" % lab
        for topo in ("N-Terminal Loss", "C-Terminal Loss"):
            g = [r for r in sub if r["topology"] == topo]
            if g:
                line += "  %s %.1f%% (n=%d)" % (topo.split("-")[0], 100.0*sum(1 for r in g if r["iso_ipTM"] >= 0.8)/len(g), len(g))
        print(line)

    json.dump(out, open(os.path.join(OUT, "sensitivity.json"), "w"), indent=2)
    print("\nwritten to out/sensitivity.json")


if __name__ == "__main__":
    main()
