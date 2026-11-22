"""Matched random-window null for the feature-deletion rate.

The feature-deletion rate counts how often a deletion window overlaps a curated
UniProt feature. With a mean deletion covering a large fraction of the chain,
the rate a randomly placed window of the same length would achieve is the
quantity it has to be compared against.
"""
import csv, json, math, os, random, re, time, urllib.parse, urllib.request
from collections import Counter

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

U = os.environ.get("AFDB_INPUT_DIR", ".")
OUT, DATA = os.path.join(BASE, "out"), os.path.join(BASE, "data")
FEATFILE = os.path.join(DATA, "uniprot_features_broad.tsv")
DRAWS = 2000
SEED = 20261119

FIELDS = ["accession", "length", "ft_region", "ft_compbias", "ft_binding", "ft_domain",
          "ft_topo_dom", "ft_transmem", "ft_repeat", "ft_motif", "ft_act_site",
          "ft_coiled", "ft_zn_fing", "ft_propep", "ft_dna_bind", "ft_intramem",
          "ft_signal", "ft_transit"]
DISORDER_COLS = {"Compositional bias", "Region"}
INHIB = re.compile(r"autoinhib|auto-inhib|inhibitor|negative regulat|repress", re.I)


def fetch_features(accs):
    if os.path.exists(FEATFILE):
        return
    base = "https://rest.uniprot.org/uniprotkb/search"
    out = open(FEATFILE, "w")
    wrote = False
    for i in range(0, len(accs), 80):
        chunk = accs[i:i + 80]
        q = " OR ".join("accession:%s" % a for a in chunk)
        url = base + "?" + urllib.parse.urlencode(
            {"query": q, "format": "tsv", "fields": ",".join(FIELDS), "size": 500})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=180) as r:
                    txt = r.read().decode("utf-8", "replace")
                break
            except Exception as exc:
                if attempt == 3:
                    raise SystemExit("FATAL: feature fetch failed: %s" % exc)
                time.sleep(3 * (attempt + 1))
        lines = txt.splitlines()
        if not wrote:
            out.write(lines[0] + "\n"); wrote = True
        for l in lines[1:]:
            out.write(l + "\n")
        out.flush()
        print("  features %d/%d" % (min(i + 80, len(accs)), len(accs)), flush=True)
        time.sleep(0.2)
    out.close()


def intervals(cell):
    """Every start..end span in a UniProt feature cell."""
    if not cell:
        return []
    return [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", cell)]


def main():
    s7 = list(csv.DictReader(open(os.path.join(U, "emergence_deletion_windows.csv"))))
    accs = sorted({r["canonical_acc"] for r in s7})
    print("emergence canonicals: %d" % len(accs))
    fetch_features(accs)

    feats, feats_dis, feats_inhib = {}, {}, {}
    with open(FEATFILE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = r["Entry"]
            allsp, dis, inh = [], [], []
            for col, cell in r.items():
                if col in ("Entry", "Length") or not cell:
                    continue
                sp = intervals(cell)
                allsp += sp
                if col in DISORDER_COLS:
                    dis += sp
                if INHIB.search(cell):
                    inh += sp
            feats[acc] = allsp; feats_dis[acc] = dis; feats_inhib[acc] = inh
    print("accessions with features parsed: %d" % len(feats))

    def hits(spans, s, e):
        return any(not (b < s or a > e) for a, b in spans)

    rnd = random.Random(SEED)
    results = {}
    for label, table in (("exact-match windows", "pub"), ("alignment windows", "aln")):
        if table == "aln":
            alnrows = {(r["gene"], r["isoform_acc"]): r for r in
                       csv.DictReader(open(os.path.join(OUT, "alignment_topology.csv")))}
        obs = obs_dis = obs_inh = 0
        exp = exp_dis = exp_inh = 0.0
        n = 0
        for r in s7:
            acc = r["canonical_acc"]
            L = int(float(r["canon_len"]))
            if table == "pub":
                m = re.findall(r"(\d+)\s*[-–]\s*(\d+)", r["deletion_window"])
                if not m:
                    continue
                s, e = int(m[0][0]), int(m[0][1])
                dl = e - s + 1
            else:
                a = alnrows.get((r["gene"], r["isoform_acc"]))
                if not a or not a["deletion_window"]:
                    continue
                ws = [tuple(map(int, w.split("-"))) for w in a["deletion_window"].split(";")]
                s, e = min(x for x, _ in ws), max(y for _, y in ws)
                dl = int(a["deletion_length"]); L = int(a["canon_len"])
            n += 1
            sp, sd, si = feats.get(acc, []), feats_dis.get(acc, []), feats_inhib.get(acc, [])
            obs += hits(sp, s, e); obs_dis += hits(sd, s, e); obs_inh += hits(si, s, e)
            if L - dl + 1 <= 1:
                exp += 1.0; exp_dis += bool(sd); exp_inh += bool(si)
                continue
            h = hd = hi = 0
            for _ in range(DRAWS):
                st = rnd.randint(1, L - dl + 1)
                en = st + dl - 1
                h += hits(sp, st, en); hd += hits(sd, st, en); hi += hits(si, st, en)
            exp += h / DRAWS; exp_dis += hd / DRAWS; exp_inh += hi / DRAWS
        results[label] = {
            "n": n,
            "observed_any_feature": 100.0 * obs / n,
            "expected_any_feature": 100.0 * exp / n,
            "observed_disorder_or_compbias": 100.0 * obs_dis / n,
            "expected_disorder_or_compbias": 100.0 * exp_dis / n,
            "observed_inhibitory": 100.0 * obs_inh / n,
            "expected_inhibitory": 100.0 * exp_inh / n,
        }
        d = results[label]
        print("\n=== %s (n=%d, %d random draws per event) ===" % (label, n, DRAWS))
        print("  any curated feature deleted      observed %5.1f%%   random null %5.1f%%   (pipeline 83.2%%)"
              % (d["observed_any_feature"], d["expected_any_feature"]))
        print("  disordered / compositional bias  observed %5.1f%%   random null %5.1f%%   (pipeline 46.6%%)"
              % (d["observed_disorder_or_compbias"], d["expected_disorder_or_compbias"]))
        print("  inhibitory or regulatory wording observed %5.1f%%   random null %5.1f%%   (pipeline 4.3%%)"
              % (d["observed_inhibitory"], d["expected_inhibitory"]))

    json.dump({"draws": DRAWS, "seed": SEED, "results": results},
              open(os.path.join(OUT, "random_window_control.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
