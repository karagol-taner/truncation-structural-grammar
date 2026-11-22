"""Stream the AFDB homodimer metadata release and keep the human rows."""
import csv, io, os, sys, time, urllib.request

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

URL = ("https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/"
       "nvda/homodimer_metadata.csv")
OUT = os.path.join(BASE, "data", "human_homodimer_metadata.csv")
LOG = os.path.join(BASE, "logs", "fetch_metadata.log")
CHUNK = 128 * 1024 * 1024
TAXID = "9606"


def log(msg):
    line = "%s  %s" % (time.strftime("%H:%M:%S"), msg)
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def get(url, start, end, attempts=5):
    for a in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, end)})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as exc:
            log("chunk %d-%d attempt %d failed: %s" % (start, end, a + 1, exc))
            time.sleep(5 * (a + 1))
    raise SystemExit("FATAL: chunk %d-%d failed after %d attempts" % (start, end, attempts))


def main():
    req = urllib.request.Request(URL, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        total = int(r.headers["Content-Length"])
    log("content-length %d bytes (%.2f GB)" % (total, total / 1e9))

    header = None
    tax_i = acc_i = None
    kept = scanned = 0
    carry = b""
    out = open(OUT, "w", newline="")
    writer = None
    pos = 0

    while pos < total:
        end = min(pos + CHUNK - 1, total - 1)
        blob = carry + get(URL, pos, end)
        pos = end + 1
        nl = blob.rfind(b"\n")
        if nl == -1:
            carry = blob
            continue
        body, carry = blob[: nl + 1], blob[nl + 1 :]

        rdr = csv.reader(io.StringIO(body.decode("utf-8", "replace")))
        for row in rdr:
            if header is None:
                header = row
                tax_i = header.index("taxId")
                acc_i = header.index("uniprotAccession")
                writer = csv.writer(out)
                writer.writerow(header)
                log("columns: %d, taxId at %d, uniprotAccession at %d"
                    % (len(header), tax_i, acc_i))
                continue
            scanned += 1
            if len(row) > tax_i and row[tax_i] == TAXID:
                writer.writerow(row)
                kept += 1
        out.flush()
        log("%.1f%%  scanned %d  kept %d" % (100.0 * pos / total, scanned, kept))

    if carry.strip():
        rdr = csv.reader(io.StringIO(carry.decode("utf-8", "replace")))
        for row in rdr:
            scanned += 1
            if len(row) > tax_i and row[tax_i] == TAXID:
                writer.writerow(row)
                kept += 1
    out.close()
    if kept == 0:
        raise SystemExit("FATAL: no human rows kept")
    log("DONE scanned %d kept %d -> %s" % (scanned, kept, OUT))


if __name__ == "__main__":
    main()
