"""Fetch GO annotations for every canonical, for a self-contained enrichment test."""
import csv, os, time, urllib.parse, urllib.request

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

BASE = "https://rest.uniprot.org/uniprotkb/search"
ACCS = os.path.join(os.environ.get("AFDB_INPUT_DIR", "."), "target_accessions.txt")
OUT = os.path.join(BASE, "data", "uniprot_go.tsv")
LOG = os.path.join(BASE, "logs", "fetch_go.log")
FIELDS = "accession,gene_primary,go_id,go_p,go_f,go_c"
BATCH = 90


def log(m):
    line = "%s  %s" % (time.strftime("%H:%M:%S"), m)
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def fetch(params, attempts=5):
    url = BASE + "?" + urllib.parse.urlencode(params)
    for a in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:
            log("attempt %d failed: %s" % (a + 1, exc))
            time.sleep(4 * (a + 1))
    raise SystemExit("FATAL: GO batch failed")


def main():
    with open(ACCS) as f:
        accs = [l.strip() for l in f if l.strip() and "-" not in l.strip()]
    log("accessions: %d" % len(accs))
    out = open(OUT, "w")
    wrote = False
    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        q = " OR ".join("accession:%s" % a for a in chunk)
        txt = fetch({"query": q, "format": "tsv", "fields": FIELDS, "size": 500})
        lines = txt.splitlines()
        if not lines:
            continue
        if not wrote:
            out.write(lines[0] + "\n"); wrote = True
        for l in lines[1:]:
            out.write(l + "\n")
        out.flush()
        if (i // BATCH) % 10 == 0:
            log("%d/%d" % (min(i + BATCH, len(accs)), len(accs)))
        time.sleep(0.2)
    out.close()
    n = sum(1 for _ in open(OUT)) - 1
    if n < 100:
        raise SystemExit("FATAL: only %d GO rows retrieved" % n)
    log("DONE %d rows -> %s" % (n, OUT))


if __name__ == "__main__":
    main()
