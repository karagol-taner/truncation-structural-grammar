"""Predicted homodimer confidence against deposited structures.

Maps every canonical accession to its PDB entries through SIFTS, reads the
deposited biological assemblies from PDBe, and records the largest number of
copies of that protein in any assembly. Two or more copies is a necessary
condition for a homodimer and is the quantity compared against the predicted
confidence tiers.

Inputs:  analysis/pairs_ipsae.csv
         SIFTS pdb_chain_uniprot.tsv.gz and the PDBe assembly endpoint
Outputs: analysis/pdb_homomer_state.csv, auxiliary/pdb_crosscheck.json
"""
import csv, gzip, json, os, sys, time
from collections import defaultdict
from math import comb

import requests

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
CACHE = os.environ.get("AFDB_PDB_CACHE",
                       os.path.join(os.environ.get("AFDB_REVISION_DIR", "."), "pdbcache"))
SIFTS = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
UNIPROT = ("https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A9606%29"
           "&fields=accession%2Cxref_pdb&format=tsv")
ASSEMBLY = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/"
BATCH = 400

os.makedirs(CACHE, exist_ok=True)


def download(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    with requests.get(url, stream=True, timeout=900) as r:
        r.raise_for_status()
        with open(path, "wb") as fh:
            for block in r.iter_content(1 << 20):
                fh.write(block)
    return path


def accession_to_entries(accessions):
    tsv = download(UNIPROT, os.path.join(CACHE, "human_xref_pdb.tsv"))
    xref = {}
    with open(tsv, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            xref[r["Entry"]] = [x for x in (r.get("PDB") or "").split(";") if x]
    out = {}
    for a in accessions:
        out[a] = xref.get(a) or xref.get(a.split("-")[0]) or []
    return out


def fetch_assemblies(pdb_ids):
    path = os.path.join(CACHE, "assemblies.jsonl")
    done = {}
    if os.path.exists(path):
        for line in open(path):
            try:
                d = json.loads(line)
                done[d["pdb_id"]] = d["assemblies"]
            except ValueError:
                pass
    todo = [p for p in pdb_ids if p not in done]
    if todo:
        s = requests.Session()
        with open(path, "a") as fh:
            for k in range(0, len(todo), BATCH):
                chunk = todo[k:k + BATCH]
                d = {}
                for attempt in range(4):
                    try:
                        r = s.post(ASSEMBLY, data=",".join(c.lower() for c in chunk),
                                   headers={"Content-Type": "application/x-www-form-urlencoded"},
                                   timeout=300)
                        r.raise_for_status()
                        d = r.json()
                        break
                    except Exception:
                        if attempt == 3:
                            raise
                        time.sleep(2 * (attempt + 1))
                for pid in chunk:
                    done[pid] = d.get(pid.lower())
                    fh.write(json.dumps({"pdb_id": pid, "assemblies": done[pid]}) + "\n")
                fh.flush()
    return done


def chain_map():
    path = download(SIFTS, os.path.join(CACHE, "pdb_chain_uniprot.tsv.gz"))
    m = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for r in csv.DictReader(fh, delimiter="\t"):
            m[(r["PDB"].lower(), r["CHAIN"])].add(r["SP_PRIMARY"])
    return m


def fisher(a, b, c, d):
    def hyp(a, b, c, d):
        return comb(a + b, a) * comb(c + d, c) / comb(a + b + c + d, a + c)
    obs, tot = hyp(a, b, c, d), 0.0
    for i in range(0, min(a + b, a + c) + 1):
        j, k, l = a + b - i, a + c - i, d - (a - i)
        if j < 0 or k < 0 or l < 0:
            continue
        p = hyp(i, j, k, l)
        if p <= obs * (1 + 1e-9):
            tot += p
    return min(1.0, tot)


def num(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError):
        return None


def main():
    pairs = list(csv.DictReader(open(os.path.join(ROOT, "analysis", "pairs_ipsae.csv"))))
    accessions = sorted({r["canonical_acc"] for r in pairs})
    acc_pdb = accession_to_entries(accessions)
    entries = sorted({p for v in acc_pdb.values() for p in v})
    asm = fetch_assemblies(entries)
    chains = chain_map()

    def state(acc):
        best = pure = seen = 0
        stem = acc.split("-")[0]
        for pid in acc_pdb.get(acc, []):
            a = asm.get(pid)
            if not a:
                continue
            seen += 1
            for assembly in a:
                mine = others = 0
                for e in assembly.get("entities", []):
                    if not str(e.get("molecule_type", "")).startswith("polypeptide"):
                        continue
                    accs = set()
                    for c in e.get("in_chains") or []:
                        accs |= {x.split("-")[0] for x in chains.get((pid.lower(), c), ())}
                    n = e.get("number_of_copies") or len(e.get("in_chains") or [])
                    if stem in accs:
                        mine += n
                    elif accs:
                        others += n
                best = max(best, mine)
                if others == 0:
                    pure = max(pure, mine)
        return best, pure, seen

    first = {}
    for r in pairs:
        first.setdefault(r["canonical_acc"], r)
    rows = []
    for a, r in first.items():
        b, p, n = state(a)
        rows.append({"canonical_acc": a, "gene": r["gene"],
                     "canon_ipTM": r["canon_ipTM"], "canon_ipSAE_min": r["canon_ipSAE_min"],
                     "n_pdb_entries": n, "max_copies": b, "max_copies_homomeric": p})
    with open(os.path.join(ROOT, "analysis", "pdb_homomer_state.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    st = {r["canonical_acc"]: r for r in rows}
    have = [r for r in rows if r["n_pdb_entries"] > 0]

    def rate(sub, key="max_copies"):
        k = sum(1 for r in sub if r[key] >= 2)
        return {"n": len(sub), "k": k,
                "pct": round(100.0 * k / len(sub), 1) if sub else None}

    tiers = {}
    for label, test in (("high", lambda v: v >= 0.8),
                        ("moderate", lambda v: 0.6 <= v < 0.8),
                        ("low", lambda v: v < 0.6)):
        sub = [r for r in have if num(r, "canon_ipTM") is not None and test(num(r, "canon_ipTM"))]
        tiers[label] = rate(sub)
        tiers[label + "_homomeric"] = rate(sub, "max_copies_homomeric")
    ips = {}
    for label, test in (("ipsae_high", lambda v: v >= 0.6), ("ipsae_low", lambda v: v < 0.6)):
        sub = [r for r in have if num(r, "canon_ipSAE_min") is not None and test(num(r, "canon_ipSAE_min"))]
        ips[label] = rate(sub)
        ips[label + "_homomeric"] = rate(sub, "max_copies_homomeric")

    em = [r for r in pairs if r["in_emergence_cohort"] in ("True", "true", "1")]
    rep = [r for r in em
           if num(r, "canon_ipSAE_min") is not None and num(r, "iso_ipSAE_min") is not None
           and num(r, "canon_ipSAE_min") < 0.6 and num(r, "iso_ipSAE_min") >= 0.6]
    repset = {id(r) for r in rep}
    non = [r for r in em if id(r) not in repset]

    def cohort(group):
        accs = {r["canonical_acc"] for r in group}
        wp = [st[a] for a in accs if st[a]["n_pdb_entries"] > 0]
        k = sum(1 for r in wp if r["max_copies"] >= 2)
        return {"events": len(group), "canonicals": len(accs), "with_structure": len(wp),
                "homomer": k, "pct": round(100.0 * k / len(wp), 1) if wp else None,
                "monomeric_genes": sorted({r["gene"] for r in wp if r["max_copies"] < 2})}

    res = {
        "accessions": len(rows),
        "accessions_with_structure": len(have),
        "pdb_entries": len(entries),
        "overall": rate(have),
        "overall_homomeric": rate(have, "max_copies_homomeric"),
        "by_ipTM_tier": tiers,
        "by_ipSAE": ips,
        "emergence_cohort": cohort(em),
        "emergence_reproducing": cohort(rep),
        "emergence_not_reproducing": cohort(non),
    }
    # the emergence cohort against other canonicals of the same low confidence,
    # stratified by how many entries each protein has, since the statistic is a
    # maximum over assemblies and rises with coverage
    em_acc = {r["canonical_acc"] for r in em}
    low_acc = {r["canonical_acc"] for r in pairs
               if num(r, "canon_ipTM") is not None and num(r, "canon_ipTM") < 0.6}

    def with_structure(accs):
        return [st[a] for a in accs if a in st and st[a]["n_pdb_entries"] > 0]

    gA, gB = with_structure(em_acc), with_structure(low_acc - em_acc)
    kA = sum(1 for r in gA if r["max_copies"] >= 2)
    kB = sum(1 for r in gB if r["max_copies"] >= 2)
    strata = []
    for lo, hi, label in ((1, 1, "1"), (2, 3, "2-3"), (4, 9, "4-9"), (10, 10 ** 9, "10+")):
        a = [r for r in gA if lo <= r["n_pdb_entries"] <= hi]
        b = [r for r in gB if lo <= r["n_pdb_entries"] <= hi]
        strata.append({"entries": label,
                       "emergence": rate(a), "other_low_confidence": rate(b)})
    res["emergence_vs_other_low"] = {
        "emergence": {"n": len(gA), "k": kA, "pct": round(100.0 * kA / len(gA), 1)},
        "other_low_confidence": {"n": len(gB), "k": kB, "pct": round(100.0 * kB / len(gB), 1)},
        "fisher": fisher(kA, len(gA) - kA, kB, len(gB) - kB),
        "by_entry_count": strata,
    }

    res["fisher_high_vs_low_ipTM"] = fisher(
        tiers["high"]["k"], tiers["high"]["n"] - tiers["high"]["k"],
        tiers["low"]["k"], tiers["low"]["n"] - tiers["low"]["k"])
    res["fisher_ipSAE"] = fisher(
        ips["ipsae_high"]["k"], ips["ipsae_high"]["n"] - ips["ipsae_high"]["k"],
        ips["ipsae_low"]["k"], ips["ipsae_low"]["n"] - ips["ipsae_low"]["k"])
    r1, r2 = res["emergence_reproducing"], res["emergence_not_reproducing"]
    res["fisher_reproducing_vs_not"] = fisher(
        r1["homomer"], r1["with_structure"] - r1["homomer"],
        r2["homomer"], r2["with_structure"] - r2["homomer"])

    with open(os.path.join(ROOT, "auxiliary", "pdb_crosscheck.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


main()
