"""The graded-response claim tested on interface size instead of on confidence.

The rheostat claim reads a graded effect off ipTM. ipTM reports
confidence in a predicted interface and not its extent, so the release's
interface residue counts are the quantity the claim needs.
"""
import csv, json, math, os
from collections import Counter

from scipy import stats

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT = os.path.join(BASE, "out")


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def med(v):
    v = sorted(x for x in v if not math.isnan(x))
    return v[len(v) // 2] if v else float("nan")


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"))))
    for r in rows:
        for k in ("canon_ipTM", "iso_ipTM", "canon_ipSAE_min", "iso_ipSAE_min",
                  "canon_iface_residues", "iso_iface_residues", "d_ipTM", "d_ipSAE",
                  "d_iface_residues", "canon_contacts", "iso_contacts"):
            r[k] = f(r[k])
        for k in ("canon_len", "iso_len", "deletion_length"):
            try: r[k] = int(float(r[k]))
            except (TypeError, ValueError): r[k] = 0
    out = {}

    real = [r for r in rows
            if r["canon_ipSAE_min"] >= 0.6 and r["canon_len"] > 0
            and r["canon_iface_residues"] > 0 and not math.isnan(r["iso_iface_residues"])]
    print("pairs whose canonical has a real interface (ipSAE_min >= 0.6): %d" % len(real))
    for r in real:
        r["frac_deleted"] = r["deletion_length"] / r["canon_len"]
        r["iface_ratio"] = r["iso_iface_residues"] / r["canon_iface_residues"]

    print("\n=== does deletion magnitude grade the interface, or only the confidence? ===")
    for xlab, xs in (("fraction of chain deleted", [r["frac_deleted"] for r in real]),
                     ("residues deleted", [float(r["deletion_length"]) for r in real])):
        for ylab, ys in (("interface residues retained (ratio)", [r["iface_ratio"] for r in real]),
                         ("change in ipTM", [r["d_ipTM"] for r in real]),
                         ("change in ipSAE_min", [r["d_ipSAE"] for r in real])):
            pts = [(a, b) for a, b in zip(xs, ys) if not (math.isnan(a) or math.isnan(b))]
            rho = stats.spearmanr([p[0] for p in pts], [p[1] for p in pts])
            print("  %-26s vs %-36s rho %+.3f  p %.2e  n %d"
                  % (xlab, ylab, rho.statistic, rho.pvalue, len(pts)))
    out["spearman_note"] = "frac_deleted vs interface ratio, ipTM change, ipSAE change"

    print("\n=== interface retained, binned by how much of the chain is removed ===")
    bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 1.01)]
    binned = []
    print("  %-14s %6s %10s %10s %10s %10s"
          % ("deleted", "n", "iface med", "iface IQR", "d_ipTM", "d_ipSAE"))
    for lo, hi in bins:
        sub = [r for r in real if lo <= r["frac_deleted"] < hi]
        if len(sub) < 5:
            continue
        ratios = sorted(r["iface_ratio"] for r in sub)
        q1, q3 = ratios[len(ratios)//4], ratios[3*len(ratios)//4]
        row = {"bin": [lo, hi], "n": len(sub),
               "iface_ratio_median": med(ratios), "iqr": [q1, q3],
               "d_ipTM_median": med([r["d_ipTM"] for r in sub]),
               "d_ipSAE_median": med([r["d_ipSAE"] for r in sub])}
        binned.append(row)
        print("  %-14s %6d %10.3f %10s %10.3f %10.3f"
              % ("%.0f-%.0f%%" % (100*lo, 100*min(hi, 1.0)), len(sub), row["iface_ratio_median"],
                 "%.2f-%.2f" % (q1, q3), row["d_ipTM_median"], row["d_ipSAE_median"]))
    out["binned_interface_retention"] = binned

    print("\n=== does the confidence change track the interface-size change? ===")
    for lab, key in (("ipTM", "d_ipTM"), ("ipSAE_min", "d_ipSAE")):
        pts = [(r[key], r["d_iface_residues"]) for r in real
               if not (math.isnan(r[key]) or math.isnan(r["d_iface_residues"]))]
        rho = stats.spearmanr([p[0] for p in pts], [p[1] for p in pts])
        out["%s_vs_iface_change" % lab] = {"rho": rho.statistic, "p": rho.pvalue, "n": len(pts)}
        print("  change in %-10s vs change in interface residues: rho %+.3f  p %.2e  n %d"
              % (lab, rho.statistic, rho.pvalue, len(pts)))

    print("\n=== the three high-tier outcome classes, on interface size ===")
    hi = [r for r in real if r["canon_ipTM"] >= 0.8]
    def cls(r):
        if r["iso_ipTM"] >= 0.8: return "Conservation"
        if r["iso_ipTM"] >= 0.6: return "Destabilization"
        return "Ablation"
    grp = {}
    for r in hi:
        grp.setdefault(cls(r), []).append(r)
    for k in ("Conservation", "Destabilization", "Ablation"):
        sub = grp.get(k, [])
        if not sub: continue
        print("  %-16s n=%3d   median interface retained %.3f   median residues deleted %4d"
              % (k, len(sub), med([r["iface_ratio"] for r in sub]),
                 med([float(r["deletion_length"]) for r in sub])))
        out.setdefault("high_tier_classes", {})[k] = {
            "n": len(sub), "iface_ratio_median": med([r["iface_ratio"] for r in sub]),
            "deleted_median": med([float(r["deletion_length"]) for r in sub])}
    vals = [[r["iface_ratio"] for r in grp.get(k, [])] for k in ("Conservation", "Destabilization", "Ablation")]
    if all(len(v) > 2 for v in vals):
        kw = stats.kruskal(*vals)
        print("  Kruskal-Wallis across the three classes on interface retained: H %.2f  p %.2e"
              % (kw.statistic, kw.pvalue))
        out["kruskal_iface_by_class"] = {"H": kw.statistic, "p": kw.pvalue}

    json.dump(out, open(os.path.join(OUT, "rheostat_interface.json"), "w"), indent=2)
    print("\nwritten to out/rheostat_interface.json")


if __name__ == "__main__":
    main()
