"""The four-class landscape on ipSAE_min, and the length-dilution model."""
import csv, json, math, os
from collections import Counter

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

OUT = os.path.join(BASE, "out")
rows = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"))))
num = ("canon_ipTM", "iso_ipTM", "canon_ipSAE_min", "iso_ipSAE_min", "d_ipTM", "d_ipSAE",
       "canon_contacts", "iso_contacts", "d_contacts", "canon_iface_residues",
       "iso_iface_residues", "d_iface_residues", "canon_pDockQ", "iso_pDockQ")
for r in rows:
    for k in num:
        try: r[k] = float(r[k])
        except (TypeError, ValueError): r[k] = float("nan")
    for k in ("canon_len", "iso_len", "deletion_length"):
        try: r[k] = int(float(r[k]))
        except (TypeError, ValueError): r[k] = 0
    r["in_emergence_cohort"] = r["in_emergence_cohort"] == "True"

def med(v):
    v = sorted(x for x in v if not math.isnan(x))
    return v[len(v) // 2] if v else float("nan")

def pearson(xs, ys):
    pts = [(a, b) for a, b in zip(xs, ys) if not (math.isnan(a) or math.isnan(b))]
    n = len(pts)
    if n < 3: return float("nan")
    mx, my = sum(p[0] for p in pts)/n, sum(p[1] for p in pts)/n
    num_ = sum((a-mx)*(b-my) for a, b in pts)
    den = math.sqrt(sum((a-mx)**2 for a, _ in pts)*sum((b-my)**2 for _, b in pts))
    return num_/den if den else float("nan")

# ---------- the four-class landscape, same rules, on each metric ----------
def landscape(cn, io, hi=0.8, lo=0.6):
    out = Counter()
    for r in rows:
        c, i = r[cn], r[io]
        if math.isnan(c) or math.isnan(i): out["unscored"] += 1; continue
        if c >= hi:
            out["Conserved" if i >= hi else ("Partial destabilization" if i >= lo else "Complete ablation")] += 1
        elif c >= lo:
            out["Mid-tier canonical"] += 1
        else:
            out["Emergence" if i >= lo else "Monomer baseline"] += 1
    return out

print("=== four-class landscape, identical rules, two metrics (n=%d) ===" % len(rows))
lt, ls = landscape("canon_ipTM", "iso_ipTM"), landscape("canon_ipSAE_min", "iso_ipSAE_min")
keys = ["Conserved", "Partial destabilization", "Complete ablation", "Mid-tier canonical",
        "Emergence", "Monomer baseline", "unscored"]
print("  %-26s %>8s %8s" % ("class", "ipTM", "ipSAE") if False else
      "  %-26s %8s %8s" % ("class", "ipTM", "ipSAE_min"))
for k in keys:
    if lt[k] or ls[k]:
        print("  %-26s %8d %8d" % (k, lt[k], ls[k]))
hi_t = lt["Conserved"] + lt["Partial destabilization"] + lt["Complete ablation"]
hi_s = ls["Conserved"] + ls["Partial destabilization"] + ls["Complete ablation"]
print("\n  high-tier canonicals: ipTM %d, ipSAE %d" % (hi_t, hi_s))
if hi_t: print("  conservation share by ipTM : %.1f%%  (pipeline 56.4%%)" % (100*lt["Conserved"]/hi_t))
if hi_s: print("  conservation share by ipSAE: %.1f%%" % (100*ls["Conserved"]/hi_s))

# ---------- does the conserved class survive the metric change? ----------
cons_t = [r for r in rows if r["canon_ipTM"] >= 0.8 and r["iso_ipTM"] >= 0.8]
both = sum(1 for r in cons_t if r["canon_ipSAE_min"] >= 0.8 and r["iso_ipSAE_min"] >= 0.8)
print("\n=== conservation class: %d by ipTM, of which %d also conserved on ipSAE (%.0f%%) ==="
      % (len(cons_t), both, 100*both/len(cons_t) if cons_t else 0))

# ---------- the length-dilution model ----------
em = [r for r in rows if r["in_emergence_cohort"]]
print("\n=== the length-dilution model on the emergence cohort (n=%d) ===" % len(em))
print("the model: predicted iso_ipTM = min(canon_ipTM * canon_len/iso_len, cap)")
for cap in (0.90, 0.95, 1.00):
    pred = [min(r["canon_ipTM"] * r["canon_len"] / r["iso_len"], cap) for r in em if r["iso_len"]]
    obs = [r["iso_ipTM"] for r in em if r["iso_len"]]
    mae = sum(abs(a-b) for a, b in zip(obs, pred))/len(obs)
    print("  cap %.2f  predicted mean %.4f  observed mean %.4f  MAE %.4f  r %+.4f"
          % (cap, sum(pred)/len(pred), sum(obs)/len(obs), mae, pearson(obs, pred)))
print("\nnull models for comparison, same cohort:")
mobs = sum(r["iso_ipTM"] for r in em)/len(em)
print("  constant at the cohort mean      predicted mean %.4f  MAE %.4f"
      % (mobs, sum(abs(r["iso_ipTM"]-mobs) for r in em)/len(em)))
print("  length ratio alone (no ipTM term): r %+.4f"
      % pearson([r["iso_ipTM"] for r in em], [r["canon_len"]/r["iso_len"] for r in em if r["iso_len"]]))
print("  fraction of chain deleted:         r %+.4f"
      % pearson([r["d_ipTM"] for r in em], [r["deletion_length"]/r["canon_len"] for r in em if r["canon_len"]]))

# ---------- what distinguishes the 42 that survive ----------
surv = [r for r in em if r["canon_ipSAE_min"] < 0.6 <= r["iso_ipSAE_min"]]
fell = [r for r in em if not (r["canon_ipSAE_min"] < 0.6 <= r["iso_ipSAE_min"])]
print("\n=== the %d events that reproduce on ipSAE against the %d that do not ===" % (len(surv), len(fell)))
for lab, grp in (("reproduce", surv), ("do not", fell)):
    print("  %-10s median d_ipTM %+.3f  d_ipSAE %+.3f  d_contacts %+7.1f  d_iface_res %+6.1f  deletion %4d aa  canon_ipTM %.2f"
          % (lab, med([r["d_ipTM"] for r in grp]), med([r["d_ipSAE"] for r in grp]),
             med([r["d_contacts"] for r in grp]), med([r["d_iface_residues"] for r in grp]),
             med([r["deletion_length"] for r in grp]), med([r["canon_ipTM"] for r in grp])))
print("  cleavage-site deletions among them: reproduce %d/%d, do not %d/%d"
      % (sum(1 for r in surv if r["cleavage_site_overlap"]), len(surv),
         sum(1 for r in fell if r["cleavage_site_overlap"]), len(fell)))
print("  topology of the reproducing set:", dict(Counter(r["topology_alignment"] for r in surv)))

with open(os.path.join(OUT, "emergence_surviving_ipsae.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(surv)
json.dump({"landscape_ipTM": dict(lt), "landscape_ipSAE": dict(ls),
           "conservation_surviving": both, "conservation_n": len(cons_t),
           "emergence_surviving": len(surv)},
          open(os.path.join(OUT, "landscape.json"), "w"), indent=2)
print("\nsurviving events written to out/emergence_surviving_ipsae.csv")
