"""Fetch UniProt sequences (with isoforms) and cleavage-site features."""
import os, sys, time, urllib.parse, urllib.request

BASE = os.environ.get("AFDB_REVISION_DIR", ".")

BASE = "https://rest.uniprot.org/uniprotkb/search"
ACCS = os.path.join(os.environ.get("AFDB_INPUT_DIR", "."), "target_accessions.txt")
OUTF = os.path.join(BASE, "data", "uniprot_features.tsv")
OUTS = os.path.join(BASE, "data", "uniprot_sequences.fasta")
LOG = os.path.join(BASE, "logs", "fetch_uniprot.log")
FIELDS = ("accession,id,protein_name,gene_primary,length,ft_transit,ft_signal,"
          "ft_propep,ft_init_met,cc_alternative_products")
BATCH = 90


def log(msg):
    line = "%s  %s" % (time.strftime("%H:%M:%S"), msg)
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
            log("attempt %d failed (%s): %s" % (a + 1, params["format"], exc))
            time.sleep(4 * (a + 1))
    raise SystemExit("FATAL: UniProt batch failed after %d attempts" % attempts)


def main():
    with open(ACCS) as f:
        accs = [l.strip() for l in f if l.strip() and "-" not in l.strip()]
    log("base accessions to fetch: %d" % len(accs))

    ft = open(OUTF, "w")
    sq = open(OUTS, "w")
    wrote_header = False
    n_seq = 0

    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        q = " OR ".join("accession:%s" % a for a in chunk)

        tsv = fetch({"query": q, "format": "tsv", "fields": FIELDS, "size": 500})
        lines = tsv.splitlines()
        if lines:
            if not wrote_header:
                ft.write(lines[0] + "\n")
                wrote_header = True
            for l in lines[1:]:
                ft.write(l + "\n")

        fa = fetch({"query": q, "format": "fasta", "includeIsoform": "true", "size": 500})
        sq.write(fa if fa.endswith("\n") else fa + "\n")
        n_seq += fa.count(">")

        ft.flush(); sq.flush()
        log("%d/%d accessions, %d sequences so far" % (min(i + BATCH, len(accs)), len(accs), n_seq))
        time.sleep(0.2)

    ft.close(); sq.close()
    if n_seq == 0:
        raise SystemExit("FATAL: no sequences retrieved")
    log("DONE %d sequences -> %s" % (n_seq, OUTS))


if __name__ == "__main__":
    main()
