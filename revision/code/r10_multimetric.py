"""Agreement across the four interface metrics the release distributes.

ipTM, ipSAE, pDockQ and LIS are on different scales, so absolute tiers cannot be
carried between them. Direction of change is scale-free, and quantile-matched
thresholds put each metric on the same footing as the ipTM tiers.
"""
import csv, json, math, os
from collections import Counter

from scipy import stats

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT = os.path.join(BASE, "out")
METRICS = [("ipTM", "canon_ipTM", "iso_ipTM"),
           ("ipSAE_min", "canon_ipSAE_min", "iso_ipSAE_min"),
           ("pDockQ", "canon_pDockQ", "iso_pDockQ"),
           ("LIS_min", "canon_LIS_min", "iso_LIS_min")]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def quantile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    i = min(int(q * (len(sorted_vals) - 1)), len(sorted_vals) - 1)
    return sorted_vals[i]


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"))))
    for r in rows:
        for _, c, i in METRICS:
            r[c], r[i] = f(r[c]), f(r[i])
        r["in_emergence_cohort"] = r["in_emergence_cohort"] == "True"
    out = {}

    # --- quantile-matched thresholds, anchored on the ipTM tiers ---
    canon_iptm = sorted(r["canon_ipTM"] for r in rows if not math.isnan(r["canon_ipTM"]))
    n = len(canon_iptm)
    q_hi = sum(1 for v in canon_iptm if v < 0.8) / n
    q_lo = sum(1 for v in canon_iptm if v < 0.6) / n
    print("ipTM tiers correspond to canonical quantiles %.4f (0.8) and %.4f (0.6)"
          % (q_hi, q_lo))
    thr = {}
    for name, c, _ in METRICS:
        vals = sorted(r[c] for r in rows if not math.isnan(r[c]))
        thr[name] = (quantile(vals, q_lo), quantile(vals, q_hi))
        print("  %-10s quantile-matched thresholds: low %.4f  high %.4f" % (name, *thr[name]))
    out["thresholds"] = {k: list(v) for k, v in thr.items()}

    def classify(cv, iv, lo, hi):
        if math.isnan(cv) or math.isnan(iv):
            return "unscored"
        if cv >= hi:
            return "Conservation" if iv >= hi else ("Destabilization" if iv >= lo else "Ablation")
        if cv >= lo:
            return "Mid-tier"
        return "Emergence" if iv >= lo else "Baseline"

    print("\n=== class counts under quantile-matched thresholds ===")
    cls = {}
    for name, c, i in METRICS:
        lo, hi = thr[name]
        cls[name] = [classify(r[c], r[i], lo, hi) for r in rows]
        print("  %-10s %s" % (name, dict(Counter(cls[name]))))
    out["class_counts"] = {k: dict(Counter(v)) for k, v in cls.items()}

    # --- where ipTM's classes go under ipSAE ---
    print("\n=== transition of ipTM classes under ipSAE, quantile-matched ===")
    order = ["Conservation", "Destabilization", "Ablation", "Emergence", "Mid-tier", "Baseline"]
    mat = {a: Counter() for a in order}
    for a, b in zip(cls["ipTM"], cls["ipSAE_min"]):
        if a in mat:
            mat[a][b] += 1
    hdr = "  %-16s" % "ipTM class" + "".join("%12s" % o[:11] for o in order) + "%8s" % "n"
    print(hdr)
    for a in order:
        tot = sum(mat[a].values())
        if not tot:
            continue
        print("  %-16s" % a + "".join("%12d" % mat[a][o] for o in order) + "%8d" % tot)
    out["transition_ipTM_to_ipSAE"] = {a: dict(mat[a]) for a in order}

    # --- direction of change, which needs no threshold at all ---
    em = [r for r in rows if r["in_emergence_cohort"]]
    print("\n=== direction of change across metrics, emergence cohort (n=%d) ===" % len(em))
    agree = Counter()
    for r in em:
        k = sum(1 for name, c, i in METRICS
                if not math.isnan(r[c]) and not math.isnan(r[i]) and r[i] > r[c])
        agree[k] += 1
    for k in sorted(agree, reverse=True):
        print("  increases on %d of 4 metrics: %3d events (%.1f%%)" % (k, agree[k], 100*agree[k]/len(em)))
    out["emergence_metric_agreement"] = dict(agree)

    for name, c, i in METRICS:
        up = sum(1 for r in em if not math.isnan(r[c]) and not math.isnan(r[i]) and r[i] > r[c])
        print("  %-10s increases in %3d of %d" % (name, up, len(em)))

    # --- do the metric deltas agree with one another across all pairs ---
    print("\n=== Spearman correlation of the score changes, all pairs ===")
    deltas = {name: [r[i] - r[c] for r in rows] for name, c, i in METRICS}
    names = [m[0] for m in METRICS]
    print("  %-12s" % "" + "".join("%12s" % n for n in names))
    corr = {}
    for a in names:
        line = "  %-12s" % a
        corr[a] = {}
        for b in names:
            pts = [(x, y) for x, y in zip(deltas[a], deltas[b])
                   if not (math.isnan(x) or math.isnan(y))]
            rho = stats.spearmanr([p[0] for p in pts], [p[1] for p in pts]).statistic
            corr[a][b] = rho
            line += "%12.3f" % rho
        print(line)
    out["delta_spearman"] = corr

    json.dump(out, open(os.path.join(OUT, "multimetric.json"), "w"), indent=2)
    print("\nwritten to out/multimetric.json")


if __name__ == "__main__":
    main()
