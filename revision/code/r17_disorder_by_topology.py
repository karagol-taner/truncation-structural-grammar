"""Disorder and compositional-bias deletion rate within each deletion topology."""
import csv, json, os, re, time, urllib.parse, urllib.request
from math import comb

BASE = os.environ.get("AFDB_REVISION_DIR", ".")
OUT = os.path.join(BASE, "out")
DISORDER_FIELDS = "accession,ft_region,ft_compbias"


def emergence_events():
    rows = list(csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv"), encoding="utf-8")))
    pairs = list(csv.DictReader(open(os.path.join(OUT, "pairs_ipsae.csv"), encoding="utf-8")))
    keep = {(r["canonical_acc"], r["isoform_acc"]) for r in pairs
            if str(r.get("in_emergence_cohort")) == "True"}
    return [r for r in rows if (r["canonical_acc"], r["isoform_acc"]) in keep]


def fetch(accs):
    out = {}
    for i in range(0, len(accs), 60):
        q = " OR ".join("accession:%s" % a for a in accs[i:i + 60])
        url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
            {"query": q, "format": "tsv", "fields": DISORDER_FIELDS})
        with urllib.request.urlopen(url, timeout=120) as fh:
            txt = fh.read().decode()
        for line in txt.splitlines()[1:]:
            c = line.split("\t")
            if c:
                out[c[0]] = " ".join(c[1:])
        time.sleep(0.3)
    return out


def fisher_two_sided(a, b, c, d):
    n1, n2, k, N = a + b, c + d, a + c, a + b + c + d
    pmf = lambda x: comb(n1, x) * comb(n2, k - x) / comb(N, k)
    p_obs = pmf(a)
    return sum(pmf(x) for x in range(max(0, k - n2), min(k, n1) + 1)
               if pmf(x) <= p_obs * (1 + 1e-9))


def main():
    ev = emergence_events()
    feats = fetch(sorted({r["canonical_acc"] for r in ev}))
    spans = lambda cell: [(int(x), int(y)) for x, y in re.findall(r"(\d+)\.\.(\d+)", cell or "")]
    tab = {}
    for r in ev:
        w = [tuple(map(int, x.split("-"))) for x in r["deletion_window"].split(";")
             if re.match(r"^\d+-\d+$", x)]
        if not w:
            continue
        sp = spans(feats.get(r["canonical_acc"], ""))
        hit = any(not (e < s2 or s > e2) for s, e in w for s2, e2 in sp)
        d = tab.setdefault(r["topology_alignment"], [0, 0])
        d[1] += 1
        d[0] += int(hit)

    n_hit, n_tot = tab.get("N-Terminal Loss", [0, 0])
    c_hit, c_tot = tab.get("C-Terminal Loss", [0, 0])
    res = {
        "by_topology": tab,
        "n_terminal_pct": round(100.0 * n_hit / n_tot, 1),
        "c_terminal_pct": round(100.0 * c_hit / c_tot, 1),
        "fisher_p": round(fisher_two_sided(n_hit, n_tot - n_hit, c_hit, c_tot - c_hit), 3),
        "all_pct": round(100.0 * sum(v[0] for v in tab.values())
                         / sum(v[1] for v in tab.values()), 1),
    }
    for k, v in res.items():
        if k != "by_topology":
            print("  %-16s %s" % (k, v))
    json.dump(res, open(os.path.join(OUT, "disorder_by_topology.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
