"""Threshold-margin, cleavage-site and length-confound controls on the emergence cohort."""
import csv, json, math, os

BASE = os.environ.get("AFDB_REVISION_DIR", ".")
OUT = os.path.join(BASE, "out")


def num(r, k):
    try:
        return float(r[k])
    except Exception:
        return None


def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i])
    rk = [0] * len(v)
    for p, i in enumerate(s):
        rk[i] = p
    return rk


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num_ = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num_ / den if den else 0.0


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"), encoding="utf-8")))
    em = [r for r in rows if str(r.get("in_emergence_cohort")) == "True"]

    def reproduces(r):
        c, i = num(r, "canon_ipSAE_min"), num(r, "iso_ipSAE_min")
        return c is not None and i is not None and c < 0.6 and i >= 0.6

    bands = []
    for lo, hi, lab in ((0.0, 0.02, "0.60-0.62"), (0.02, 0.05, "0.62-0.65"),
                        (0.05, 0.10, "0.65-0.70"), (0.10, 9.0, ">0.70")):
        g = [r for r in em if num(r, "iso_ipTM") is not None
             and lo <= num(r, "iso_ipTM") - 0.6 < hi]
        rep = [r for r in g if reproduces(r)]
        bands.append({"band": lab, "n": len(g), "reproduce": len(rep),
                      "pct": round(100.0 * len(rep) / len(g), 1) if g else None})

    cl = [r for r in em if str(r.get("cleavage_site_overlap", "")).strip().lower()
          not in ("", "none", "nan", "false", "0")]
    ncl = [r for r in em if r not in cl]
    cleav = {lab: {"n": len(g), "reproduce": sum(reproduces(r) for r in g),
                   "pct": round(100.0 * sum(reproduces(r) for r in g) / len(g), 1) if g else None}
             for lab, g in (("overlaps_cleavage_site", cl), ("no_overlap", ncl))}

    X, Y, Z = [], [], []
    for r in rows:
        if not num(r, "canon_iface_residues"):
            continue
        cl_, il_ = num(r, "canon_len"), num(r, "iso_len")
        if not cl_ or il_ is None:
            continue
        X.append((cl_ - il_) / cl_)
        Y.append(num(r, "d_ipTM") or 0.0)
        Z.append(num(r, "d_iface_residues") or 0.0)
    rx, ry, rz = rank(X), rank(Y), rank(Z)
    rxy, rxz, ryz = pearson(rx, ry), pearson(rx, rz), pearson(ry, rz)
    partial = (rxy - rxz * ryz) / math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))

    res = {"threshold_margin": bands, "cleavage_site": cleav,
           "length_confound": {"n_pairs": len(X),
                               "rho_deletion_vs_d_ipTM": round(rxy, 3),
                               "rho_deletion_vs_d_iface": round(rxz, 3),
                               "rho_iface_vs_d_ipTM": round(ryz, 3),
                               "partial_deletion_vs_d_ipTM": round(partial, 3)}}
    json.dump(res, open(os.path.join(OUT, "emergence_controls.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
