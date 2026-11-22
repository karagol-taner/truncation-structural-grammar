"""Fetch AFDB pLDDT disorder fractions for every canonical, for a matched control."""
import csv, json, os, time, urllib.request

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

ACCS = os.path.join(os.environ.get("AFDB_INPUT_DIR", "."), "target_accessions.txt")
OUT = os.path.join(BASE, "data", "plddt_disorder.csv")
LOG = os.path.join(BASE, "logs", "fetch_plddt.log")
FIELDS = ("uniprotAccession", "globalMetricValue", "fractionPlddtVeryLow",
          "fractionPlddtLow", "fractionPlddtConfident", "fractionPlddtVeryHigh",
          "sequenceEnd", "latestVersion")


def log(m):
    line = "%s  %s" % (time.strftime("%H:%M:%S"), m)
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def main():
    with open(ACCS) as f:
        accs = [l.strip() for l in f if l.strip() and "-" not in l.strip()]
    done = set()
    if os.path.exists(OUT):
        with open(OUT) as f:
            done = {r["uniprotAccession"] for r in csv.DictReader(f)}
    todo = [a for a in accs if a not in done]
    log("accessions %d, already fetched %d, to fetch %d" % (len(accs), len(done), len(todo)))

    new = not done
    fh = open(OUT, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=list(FIELDS) + ["status"])
    if new:
        w.writeheader()

    ok = miss = 0
    for n, a in enumerate(todo, 1):
        url = "https://alphafold.ebi.ac.uk/api/prediction/%s" % a
        rec = {k: "" for k in FIELDS}
        rec["uniprotAccession"] = a
        rec["status"] = "error"
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    d = json.loads(r.read().decode())
                e = d[0] if isinstance(d, list) and d else d
                if e:
                    for k in FIELDS:
                        rec[k] = e.get(k, "")
                    rec["uniprotAccession"] = a
                    rec["status"] = "ok"
                    ok += 1
                else:
                    rec["status"] = "empty"
                    miss += 1
                break
            except Exception as exc:
                if attempt == 3:
                    rec["status"] = "error:%s" % str(exc)[:40]
                    miss += 1
                else:
                    time.sleep(2 * (attempt + 1))
        w.writerow(rec)
        if n % 250 == 0:
            fh.flush()
            log("%d/%d  ok %d  missing %d" % (n, len(todo), ok, miss))
    fh.close()
    log("DONE ok %d missing %d" % (ok, miss))
    if ok == 0 and todo:
        raise SystemExit("FATAL: no pLDDT records retrieved")


if __name__ == "__main__":
    main()
