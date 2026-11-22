"""The disorder comparison against an ipTM-matched control.

The unmatched comparison sets the emergence cohort, drawn from canonicals below
the tier boundary, against the ablation cohort, drawn from canonicals at 0.8 and
above. Disorder depresses ipTM, so that contrast is partly the selection rule.
The comparison that bears on the claim holds the canonical ipTM stratum fixed.
"""
import csv, json, math, os
from collections import Counter

from scipy import stats

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

U = os.environ.get("AFDB_INPUT_DIR", ".")
OUT, DATA = os.path.join(BASE, "out"), os.path.join(BASE, "data")


def cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float("nan")
    va, vb = stats.tvar(a), stats.tvar(b)
    s = math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    return (sum(a) / na - sum(b) / nb) / s if s else float("nan")


def report(name, a, b, la, lb):
    if len(a) < 2 or len(b) < 2:
        print("  %-46s insufficient n (%d vs %d)" % (name, len(a), len(b)))
        return None
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cohens_d(a, b)
    print("  %-46s %s %.4f (n=%d)  %s %.4f (n=%d)   d=%+.2f  p=%.2e"
          % (name, la, sum(a) / len(a), len(a), lb, sum(b) / len(b), len(b), d, u.pvalue))
    return {"name": name, "mean_a": sum(a) / len(a), "n_a": len(a),
            "mean_b": sum(b) / len(b), "n_b": len(b), "cohens_d": d, "p": u.pvalue}


def main():
    dis = {}
    with open(os.path.join(DATA, "plddt_disorder.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["status"] == "ok" and r["fractionPlddtVeryLow"] not in ("", None):
                dis[r["uniprotAccession"]] = float(r["fractionPlddtVeryLow"])
    print("canonicals with a pLDDT disorder fraction: %d" % len(dis))

    pairs = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"))))
    for r in pairs:
        for k in ("canon_ipTM", "iso_ipTM", "canon_ipSAE_min", "iso_ipSAE_min"):
            try: r[k] = float(r[k])
            except (TypeError, ValueError): r[k] = float("nan")
        r["in_emergence_cohort"] = r["in_emergence_cohort"] == "True"
        r["in_ablation_cohort"] = r["in_ablation_cohort"] == "True"
        r["disorder"] = dis.get(r["canonical_acc"])

    have = [r for r in pairs if r["disorder"] is not None and not math.isnan(r["canon_ipTM"])]
    print("pairs with disorder available: %d of %d" % (len(have), len(pairs)))

    em = [r for r in have if r["in_emergence_cohort"]]
    ab = [r for r in have if r["in_ablation_cohort"]]
    low = [r for r in have if r["canon_ipTM"] < 0.6]
    low_non_em = [r for r in low if not r["in_emergence_cohort"]]

    out = {}
    print("\n=== the unmatched comparison, reproduced ===")
    out["unmatched"] = report("emergence vs complete-ablation cohort",
                              [r["disorder"] for r in em], [r["disorder"] for r in ab],
                              "emergence", "ablation")
    print("  pipeline figure for this contrast: Cohen's d = 1.39, p = 4.3e-17")

    print("\n=== the matched comparison: both groups drawn from canonicals below 0.6 ===")
    out["matched"] = report("emergence vs other sub-0.6 canonicals",
                            [r["disorder"] for r in em], [r["disorder"] for r in low_non_em],
                            "emergence", "control")

    print("\n=== within narrow canonical ipTM bins ===")
    bins = [(0.0, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6)]
    binned = []
    for lo, hi in bins:
        a = [r["disorder"] for r in em if lo <= r["canon_ipTM"] < hi]
        b = [r["disorder"] for r in low_non_em if lo <= r["canon_ipTM"] < hi]
        res = report("ipTM %.1f to %.1f" % (lo, hi), a, b, "emergence", "control")
        if res: binned.append({"bin": [lo, hi], **res})
    out["binned"] = binned

    print("\n=== disorder against the canonical ipTM itself, which is the confound ===")
    x = [r["canon_ipTM"] for r in have]
    y = [r["disorder"] for r in have]
    rho = stats.spearmanr(x, y)
    print("  Spearman rho = %+.4f  p = %.2e  (n=%d)" % (rho.statistic, rho.pvalue, len(x)))
    out["disorder_vs_ipTM_spearman"] = {"rho": rho.statistic, "p": rho.pvalue, "n": len(x)}

    print("\n=== do the events that reproduce on ipSAE differ in disorder? ===")
    surv = [r for r in em if r["canon_ipSAE_min"] < 0.6 <= r["iso_ipSAE_min"]]
    fell = [r for r in em if r not in surv]
    out["surviving_vs_not"] = report("ipSAE-reproducing vs not, within emergence",
                                     [r["disorder"] for r in surv], [r["disorder"] for r in fell],
                                     "reproduce", "do not")

    print("\n=== the AIUPred peak score against sequence length, in the pipeline table ===")
    s6 = list(csv.DictReader(open(os.path.join(U, "disorder_scores.csv"))))
    xs, ys, coh = [], [], []
    for r in s6:
        try:
            xs.append(float(r["Sequence_Length"])); ys.append(float(r["Peak_AIUPred_Binding_Score"]))
            coh.append(r.get("Cohort", ""))
        except (TypeError, ValueError):
            pass
    if len(xs) > 3:
        rr = stats.spearmanr(xs, ys)
        print("  peak score vs sequence length: Spearman rho = %+.4f  p = %.2e  (n=%d)"
              % (rr.statistic, rr.pvalue, len(xs)))
        for c in sorted(set(coh)):
            lens = [x for x, k in zip(xs, coh) if k == c]
            if lens:
                print("    cohort %-28s n=%3d  mean length %6.1f" % (c, len(lens), sum(lens)/len(lens)))
        out["aiupred_peak_vs_length"] = {"rho": rr.statistic, "p": rr.pvalue, "n": len(xs)}

    json.dump(out, open(os.path.join(OUT, "disorder_control.json"), "w"), indent=2, default=str)


if __name__ == "__main__":
    main()
