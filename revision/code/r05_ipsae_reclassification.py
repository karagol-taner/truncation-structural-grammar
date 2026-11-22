"""Rebuild the truncation landscape on the metric the release itself uses.

The release defines its high-confidence set and its confidence tiers on
ipSAE_min, not on ipTM, and distributes interface-contact counts alongside.
This reclassifies every pair on ipSAE_min and reports whether the predicted
interface itself changes.
"""
import csv, json, os, math
from collections import Counter, defaultdict

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

U = os.environ.get("AFDB_INPUT_DIR", ".")
DATA, OUT = os.path.join(BASE, "data"), os.path.join(BASE, "out")
TIERS = (0.6, 0.7, 0.8)


def f(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def tier(v):
    if math.isnan(v):
        return "unscored"
    if v >= 0.8:
        return "very high"
    if v >= 0.7:
        return "confident"
    if v >= 0.6:
        return "low confidence"
    return "below release threshold"


def pearson(xs, ys):
    pts = [(a, b) for a, b in zip(xs, ys) if not (math.isnan(a) or math.isnan(b))]
    if len(pts) < 3:
        return float("nan")
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    num = sum((p[0] - mx) * (p[1] - my) for p in pts)
    den = math.sqrt(sum((p[0] - mx) ** 2 for p in pts) * sum((p[1] - my) ** 2 for p in pts))
    return num / den if den else float("nan")


def main():
    meta = {}
    with open(os.path.join(DATA, "human_homodimer_metadata.csv")) as fh:
        for r in csv.DictReader(fh):
            meta[r["uniprotAccession"]] = r

    aln = {(r["gene"], r["isoform_acc"]): r
           for r in csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv")))}
    pairs = list(csv.DictReader(open(os.path.join(U, "pairs_master.csv"))))
    s5 = {(r["gene"], r["isoform_acc"]) for r in csv.DictReader(open(os.path.join(U, "emergence_cohort.csv")))}
    s3 = {(r["gene"], r["isoform_acc"]) for r in csv.DictReader(open(os.path.join(U, "ablation_cohort.csv")))}

    rows = []
    for p in pairs:
        c, i = meta.get(p["canonical_acc"]), meta.get(p["isoform_acc"])
        if not c or not i:
            continue
        key = (p["gene"], p["isoform_acc"])
        a = aln.get(key)

        def pack(m, pre):
            sae = min(f(m["ipSAE_AB"]), f(m["ipSAE_BA"]))
            contacts = f(m["numberOfInteractions"])
            iface = f(m["dist1_AB"]) + f(m["dist2_AB"])
            return {
                pre + "ipTM": f(m["ipTM"]), pre + "ipSAE_min": sae,
                pre + "ipSAE_AB": f(m["ipSAE_AB"]), pre + "ipSAE_BA": f(m["ipSAE_BA"]),
                pre + "pDockQ": f(m["pDockQ"]),
                pre + "LIS_min": min(f(m["LIS_AB"]), f(m["LIS_BA"])),
                pre + "contacts": contacts, pre + "iface_residues": iface,
                pre + "clash_backbone": f(m["N_clash_backbone"]),
                pre + "tier_ipSAE": tier(sae), pre + "tier_ipTM": tier(f(m["ipTM"])),
            }

        row = {"gene": p["gene"], "canonical_acc": p["canonical_acc"],
               "isoform_acc": p["isoform_acc"],
               "pipeline_outcome": p["structural_outcome"],
               "in_emergence_cohort": key in s5, "in_ablation_cohort": key in s3,
               "canon_len": a["canon_len"] if a else "", "iso_len": a["iso_len"] if a else "",
               "deletion_length": a["deletion_length"] if a else "",
               "topology_alignment": a["topology_alignment"] if a else "",
               "cleavage_site_overlap": a["cleavage_site_overlap"] if a else ""}
        row.update(pack(c, "canon_"))
        row.update(pack(i, "iso_"))
        row["d_ipTM"] = row["iso_ipTM"] - row["canon_ipTM"]
        row["d_ipSAE"] = row["iso_ipSAE_min"] - row["canon_ipSAE_min"]
        row["d_contacts"] = row["iso_contacts"] - row["canon_contacts"]
        row["d_iface_residues"] = row["iso_iface_residues"] - row["canon_iface_residues"]
        row["d_len"] = f(row["iso_len"]) - f(row["canon_len"])
        rows.append(row)

    with open(os.path.join(OUT, "pairs_ipsae.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    print("pairs with metadata on both members: %d of %d" % (len(rows), len(pairs)))
    print("\nipSAE_AB vs ipSAE_BA on canonicals, Pearson r = %.4f (homodimer, expected near 1)"
          % pearson([r["canon_ipSAE_AB"] for r in rows], [r["canon_ipSAE_BA"] for r in rows]))

    print("\n=== tier distribution of canonicals ===")
    print("  by ipTM     :", dict(Counter(r["canon_tier_ipTM"] for r in rows)))
    print("  by ipSAE_min:", dict(Counter(r["canon_tier_ipSAE"] for r in rows)))

    print("\n=== length confounding: correlation of the score change with the length change ===")
    print("  d_ipTM  vs d_len : r = %+.4f" % pearson([r["d_ipTM"] for r in rows], [r["d_len"] for r in rows]))
    print("  d_ipSAE vs d_len : r = %+.4f" % pearson([r["d_ipSAE"] for r in rows], [r["d_len"] for r in rows]))
    print("  d_contacts vs d_len: r = %+.4f" % pearson([r["d_contacts"] for r in rows], [r["d_len"] for r in rows]))

    em = [r for r in rows if r["in_emergence_cohort"]]
    ab = [r for r in rows if r["in_ablation_cohort"]]
    print("\n=== the emergence cohort on ipSAE and on interface contacts (n=%d) ===" % len(em))
    surv = [r for r in em if r["canon_ipSAE_min"] < 0.6 <= r["iso_ipSAE_min"]]
    print("  emergence reproduced on ipSAE_min (canonical <0.6, isoform >=0.6): %d" % len(surv))
    print("  isoform ipSAE_min exceeds canonical at all:                        %d" % sum(1 for r in em if r["d_ipSAE"] > 0))
    print("  interface contacts increase:                                       %d" % sum(1 for r in em if r["d_contacts"] > 0))
    print("  interface residue count increases:                                 %d" % sum(1 for r in em if r["d_iface_residues"] > 0))
    print("  both ipSAE and contacts increase:                                  %d"
          % sum(1 for r in em if r["d_ipSAE"] > 0 and r["d_contacts"] > 0))
    def med(v):
        v = sorted(x for x in v if not math.isnan(x))
        return v[len(v) // 2] if v else float("nan")
    print("  median d_ipTM %+.3f   median d_ipSAE %+.3f   median d_contacts %+.1f   median d_iface_res %+.1f"
          % (med([r["d_ipTM"] for r in em]), med([r["d_ipSAE"] for r in em]),
             med([r["d_contacts"] for r in em]), med([r["d_iface_residues"] for r in em])))

    print("\n=== the complete-ablation cohort for comparison (n=%d) ===" % len(ab))
    print("  isoform ipSAE_min falls below canonical: %d" % sum(1 for r in ab if r["d_ipSAE"] < 0))
    print("  interface contacts decrease:             %d" % sum(1 for r in ab if r["d_contacts"] < 0))
    print("  median d_ipTM %+.3f   median d_ipSAE %+.3f   median d_contacts %+.1f"
          % (med([r["d_ipTM"] for r in ab]), med([r["d_ipSAE"] for r in ab]),
             med([r["d_contacts"] for r in ab])))

    summary = {
        "pairs": len(rows),
        "emergence_n": len(em),
        "emergence_reproduced_on_ipSAE": len(surv),
        "emergence_contacts_increase": sum(1 for r in em if r["d_contacts"] > 0),
        "emergence_ipSAE_and_contacts_increase": sum(1 for r in em if r["d_ipSAE"] > 0 and r["d_contacts"] > 0),
        "r_dipTM_dlen": pearson([r["d_ipTM"] for r in rows], [r["d_len"] for r in rows]),
        "r_dipSAE_dlen": pearson([r["d_ipSAE"] for r in rows], [r["d_len"] for r in rows]),
        "canon_tier_ipTM": dict(Counter(r["canon_tier_ipTM"] for r in rows)),
        "canon_tier_ipSAE": dict(Counter(r["canon_tier_ipSAE"] for r in rows)),
    }
    json.dump(summary, open(os.path.join(OUT, "ipsae_summary.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
