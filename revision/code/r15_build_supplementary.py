"""Rebuild the supplementary workbook from the revision outputs."""
import csv, json, os, sys
from collections import Counter, OrderedDict

import openpyxl
from openpyxl.styles import Font

OUT = sys.argv[1]
OLD = sys.argv[2]
DST = sys.argv[3]
FIGDIR = sys.argv[4] if len(sys.argv) > 4 else "figures"
MASTER = sys.argv[5] if len(sys.argv) > 5 else os.path.join("revision", "input", "pairs_master.csv")
FEATURES = sys.argv[6] if len(sys.argv) > 6 else os.path.join("auxiliary", "Domain_Annotations_v2_per_event.csv")
FONT = Font(name="Calibri", size=11)




def rd(name):
    with open(os.path.join(OUT, name)) as fh:
        return list(csv.DictReader(fh))


def js(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else None


def sheet(wb, name, header, rows, text_cols=()):
    ws = wb.create_sheet(name)
    ws.append(header)
    for c in ws[1]:
        c.font = FONT
    ti = [header.index(c) for c in text_cols if c in header]
    for r in rows:
        ws.append(r)
        for i, c in enumerate(ws[ws.max_row]):
            c.font = FONT
            if i in ti:
                c.number_format = "@"
    return ws


def main():
    pairs = rd("pairs_ipsae.csv")
    aln = {(r["gene"], r["isoform_acc"]): r for r in rd("alignment_topology.csv")}
    plddt = {r["uniprotAccession"]: r["fractionPlddtVeryLow"]
             for r in rd("plddt_disorder.csv") if r["status"] == "ok"}
    dis, rwc = js("disorder_control.json"), js("random_window_control.json")
    mm, rh = js("multimetric.json"), js("rheostat_interface.json")
    sse, sens = js("shared_splice_events.json"), js("sensitivity.json")
    arb = js("topology_arbitration.json")
    master = {}
    mpath = os.path.join(os.path.dirname(OUT), "..", "GitHub", "revision", "input", "pairs_master.csv")
    for cand in (os.path.join(OUT, "pairs_master.csv"), MASTER, mpath):
        if cand and os.path.exists(cand):
            master = {(r["gene"], r["isoform_acc"]): r for r in
                      csv.DictReader(open(cand))}
            break
    feat = {}
    if FEATURES and os.path.exists(FEATURES):
        feat = {(r["gene"], r["isoform_acc"]): r for r in csv.DictReader(open(FEATURES))}

    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    # ---- Table S1: every pair, both metrics, realigned topology ----
    hdr = ["gene", "canonical_acc", "isoform_acc", "canon_len", "iso_len", "deletion_length",
           "deletion_window", "n_deleted_blocks", "junction_substitutions",
           "is_exact_subsequence", "topology_alignment", "cleavage_site_overlap",
           "canon_ipTM", "iso_ipTM", "canon_ipSAE_min", "iso_ipSAE_min",
           "canon_pDockQ", "iso_pDockQ", "canon_LIS_min", "iso_LIS_min",
           "canon_iface_residues", "iso_iface_residues", "d_iface_residues",
           "canon_disorder_pLDDT50", "pipeline_outcome",
           "reference_is_uniprot_displayed_canonical"]
    rows = []
    for p in pairs:
        a = aln.get((p["gene"], p["isoform_acc"]), {})
        rows.append([p["gene"], p["canonical_acc"], p["isoform_acc"],
                     a.get("canon_len", ""), a.get("iso_len", ""), a.get("deletion_length", ""),
                     a.get("deletion_window", ""), a.get("n_deleted_blocks", ""),
                     a.get("junction_substitutions", ""), a.get("is_exact_subsequence", ""),
                     a.get("topology_alignment", ""), a.get("cleavage_site_overlap", ""),
                     p["canon_ipTM"], p["iso_ipTM"], p["canon_ipSAE_min"], p["iso_ipSAE_min"],
                     p["canon_pDockQ"], p["iso_pDockQ"], p["canon_LIS_min"], p["iso_LIS_min"],
                     p["canon_iface_residues"], p["iso_iface_residues"], p["d_iface_residues"],
                     plddt.get(p["canonical_acc"], ""), p["pipeline_outcome"],
                     master.get((p["gene"], p["isoform_acc"]), {})
                           .get("reference_is_uniprot_displayed_canonical", "")])
    sheet(wb, "Table S1", hdr, rows, text_cols=("deletion_window", "cleavage_site_overlap"))

    # ---- Table S2: topology against outcome, realigned ----
    tab = Counter()
    for p in pairs:
        a = aln.get((p["gene"], p["isoform_acc"]))
        if a:
            tab[(a["topology_alignment"], p["pipeline_outcome"])] += 1
    topo_tot = Counter()
    for (t, _), n in tab.items():
        topo_tot[t] += n
    rows = [[t, o, n, topo_tot[t], round(100.0 * n / topo_tot[t], 2)]
            for (t, o), n in sorted(tab.items())]
    sheet(wb, "Table S2", ["truncation_topology", "pipeline_outcome", "count",
                           "topology_total", "percent_of_topology"], rows)

    # ---- S3, S4, S5: the three cohorts ----
    def cohort_rows(sel):
        out = []
        for p in pairs:
            if not sel(p):
                continue
            a = aln.get((p["gene"], p["isoform_acc"]), {})
            out.append([p["gene"], p["canonical_acc"], p["isoform_acc"],
                        a.get("canon_len", ""), a.get("iso_len", ""),
                        a.get("deletion_length", ""), a.get("deletion_window", ""),
                        a.get("topology_alignment", ""), a.get("cleavage_site_overlap", ""),
                        p["canon_ipTM"], p["iso_ipTM"], p["d_ipTM"],
                        p["canon_ipSAE_min"], p["iso_ipSAE_min"], p["d_ipSAE"],
                        p["canon_iface_residues"], p["iso_iface_residues"],
                        p["d_iface_residues"], plddt.get(p["canonical_acc"], "")])
        return out

    ch = ["gene", "canonical_acc", "isoform_acc", "canon_len", "iso_len", "deletion_length",
          "deletion_window", "topology_alignment", "cleavage_site_overlap",
          "canon_ipTM", "iso_ipTM", "d_ipTM", "canon_ipSAE_min", "iso_ipSAE_min", "d_ipSAE",
          "canon_iface_residues", "iso_iface_residues", "d_iface_residues",
          "canon_disorder_pLDDT50"]
    sheet(wb, "Table S3", ch, cohort_rows(lambda p: p["in_ablation_cohort"] == "True"),
          text_cols=("deletion_window", "cleavage_site_overlap"))

    def is_destab(p):
        return float(p["canon_ipTM"]) >= 0.8 and 0.6 <= float(p["iso_ipTM"]) < 0.8
    sheet(wb, "Table S4", ch, cohort_rows(is_destab),
          text_cols=("deletion_window", "cleavage_site_overlap"))

    em = cohort_rows(lambda p: p["in_emergence_cohort"] == "True")
    emk = [(p["gene"], p["isoform_acc"]) for p in pairs if p["in_emergence_cohort"] == "True"]
    repro = {}
    for p in pairs:
        if p["in_emergence_cohort"] == "True":
            repro[(p["gene"], p["isoform_acc"])] = (
                float(p["canon_ipSAE_min"]) < 0.6 <= float(p["iso_ipSAE_min"]))
    for row, key in zip(em, emk):
        row.append("TRUE" if repro[key] else "FALSE")
    sheet(wb, "Table S5", ch + ["reproduces_on_ipSAE"], em,
          text_cols=("deletion_window", "cleavage_site_overlap"))

    # ---- Table S6: disorder, all canonicals with cohort labels ----
    lab = {}
    for p in pairs:
        a = p["canonical_acc"]
        if p["in_emergence_cohort"] == "True":
            lab[a] = "emergence"
        elif p["in_ablation_cohort"] == "True":
            lab.setdefault(a, "complete ablation")
        else:
            lab.setdefault(a, "other")
    ciptm = {}
    for p in pairs:
        ciptm.setdefault(p["canonical_acc"], p["canon_ipTM"])
    rows = [[a, lab.get(a, "other"), ciptm.get(a, ""), v] for a, v in sorted(plddt.items())]
    sheet(wb, "Table S6", ["canonical_acc", "cohort", "canon_ipTM",
                           "disorder_fraction_pLDDT50"], rows)

    # ---- Table S7: feature mining on realigned windows ----
    rows = []
    for p in pairs:
        if p["in_emergence_cohort"] != "True":
            continue
        a = aln.get((p["gene"], p["isoform_acc"]), {})
        g = feat.get((p["gene"], p["isoform_acc"]), {})
        rows.append([p["gene"], p["canonical_acc"], p["isoform_acc"],
                     a.get("canon_len", ""), a.get("iso_len", ""),
                     a.get("topology_alignment", ""), a.get("deletion_window", ""),
                     a.get("deletion_length", ""), a.get("cleavage_site_overlap", ""),
                     g.get("n_features_overlapping", ""), g.get("features_overlapping", ""),
                     g.get("any_inhibitory_lexical", ""), g.get("any_disorder_lexical", "")])
    sheet(wb, "Table S7", ["gene", "canonical_acc", "isoform_acc", "canon_len", "iso_len",
                           "topology", "deletion_window", "deletion_length",
                           "cleavage_site_overlap", "n_features_overlapping",
                           "features_overlapping", "any_inhibitory_lexical",
                           "any_disorder_lexical"], rows,
          text_cols=("deletion_window", "cleavage_site_overlap", "features_overlapping"))

    # ---- Table S8 carried over unchanged ----
    old = openpyxl.load_workbook(OLD, data_only=True)
    if "Table S8" in old.sheetnames:
        src = old["Table S8"]
        ws = wb.create_sheet("Table S8")
        for row in src.iter_rows(values_only=True):
            ws.append(row)
            for c in ws[ws.max_row]:
                c.font = FONT

    # ---- Table S9: enrichment, from the g:Profiler exports ----
    GDIR = os.path.join(OUT, "gprofiler")
    EXPORTS = [
        ("emergence_188", "custom, background including unannotated genes",
         "result_emergence_188_unannotated_background.csv"),
        ("emergence_188", "custom over annotated genes", "result_emergence_188.csv"),
        ("complete_ablation_55", "custom over annotated genes", "result_complete_ablation.csv"),
    ]
    rows = []
    for cohort, scope, fname in EXPORTS:
        fp = os.path.join(GDIR, fname)
        if not os.path.exists(fp):
            continue
        for r in csv.DictReader(open(fp, encoding="utf-8-sig")):
            ts, qs = int(r["term_size"]), int(r["query_size"])
            isz = int(r["intersection_size"])
            rows.append([cohort, scope, r["source"], r["term_id"], r["term_name"],
                         float(r["adjusted_p_value"]), ts, qs, isz,
                         int(r["effective_domain_size"]),
                         round(isz / qs, 6), round(isz / ts, 6), r["intersections"]])
    sheet(wb, "Table S9",
          ["cohort", "domain_scope", "source", "term_id", "term_name", "p_adjusted",
           "term_size", "query_size", "intersection_size", "effective_domain_size",
           "precision", "recall", "intersections"], rows)

    # ---- Table S10: metric comparison ----
    rows = []
    if mm:
        for metric, counts in mm["class_counts"].items():
            for cls, n in sorted(counts.items()):
                rows.append(["class count under quantile-matched thresholds", metric, cls, n])
        for a, inner in mm["delta_spearman"].items():
            for b, v in inner.items():
                rows.append(["Spearman of score change", a, b, round(v, 4)])
        for a, inner in mm["transition_ipTM_to_ipSAE"].items():
            for b, n in inner.items():
                rows.append(["ipTM class to ipSAE class", a, b, n])
        for k, v in sorted(mm["emergence_metric_agreement"].items()):
            rows.append(["emergence events increasing on N of 4 metrics", k, "", v])
    sheet(wb, "Table S10", ["quantity", "a", "b", "value"], rows)

    # ---- Table S11: controls ----
    rows = []
    if rwc:
        for wset, d in rwc["results"].items():
            for k, v in d.items():
                rows.append(["random-window control", wset, k,
                             round(v, 3) if isinstance(v, float) else v])
    if dis:
        for key in ("unmatched", "matched", "surviving_vs_not"):
            if dis.get(key):
                for k, v in dis[key].items():
                    rows.append(["disorder contrast", key, k,
                                 round(v, 5) if isinstance(v, float) else v])
        for b in dis.get("binned", []):
            rows.append(["disorder within canonical ipTM stratum", b["name"],
                         "cohens_d", round(b["cohens_d"], 3)])
            rows.append(["disorder within canonical ipTM stratum", b["name"],
                         "p", b["p"]])
    if rh:
        for b in rh.get("binned_interface_retention", []):
            rows.append(["interface retained by deletion fraction",
                         "%.0f-%.0f%%" % (100 * b["bin"][0], 100 * min(b["bin"][1], 1.0)),
                         "median_ratio", round(b["iface_ratio_median"], 4)])
    if sens:
        for k, v in sens.items():
            if isinstance(v, list):
                rows.append(["bootstrap 95% interval", k, "", "%.1f to %.1f" % tuple(v)])
    sheet(wb, "Table S11", ["quantity", "a", "b", "value"], rows)

    # ---- Table S12: shared splicing events and topology arbitration ----
    rows = []
    if sse:
        for label, d in sse.items():
            rows.append(["cohort", label, "events", d["events"]])
            rows.append(["cohort", label, "distinct_genes", d["distinct_genes"]])
            for thr, t in sorted(d["by_threshold"].items()):
                rows.append(["independent events at deleted-segment identity", label, thr,
                             t["independent_events"]])
                for c in t["clusters"]:
                    rows.append(["cluster", label, thr, ", ".join(c["genes"])])
    if arb:
        for tbl, counts in arb["topology_counts"].items():
            for k, v in sorted(counts.items()):
                rows.append(["emergence topology by assignment", tbl, k, v])
    sheet(wb, "Table S12", ["quantity", "a", "b", "value"], rows)

    # ---- Table S13: deposited-structure state ----
    pdbst = []
    pdbp = os.path.join(OUT, "pdb_homomer_state.csv")
    if os.path.exists(pdbp):
        for r in sorted(csv.DictReader(open(pdbp)), key=lambda x: (x["gene"], x["canonical_acc"])):
            if int(r["n_pdb_entries"]) > 0:
                pdbst.append([r["gene"], r["canonical_acc"], r["canon_ipTM"], r["canon_ipSAE_min"],
                              int(r["n_pdb_entries"]), int(r["max_copies"]),
                              int(r["max_copies_homomeric"])])
    if pdbst:
        sheet(wb, "Table S13",
              ["gene", "canonical_acc", "canon_ipTM", "canon_ipSAE_min", "n_pdb_entries",
               "max_copies_any_assembly", "max_copies_homomeric_assembly"], pdbst)

    # ---- Title page ----
    ws = wb.create_sheet("Title Page", 0)
    toc = OrderedDict([
        ("Table S1", "Main landscape dataset: every pair on four interface metrics, with realigned deletion windows"),
        ("Table S2", "Aggregated topology against outcome class"),
        ("Table S3", "Complete-ablation cohort"),
        ("Table S4", "Partial-destabilization cohort"),
        ("Table S5", "Emergence cohort, with reproduction on ipSAE flagged"),
        ("Table S6", "Canonical pLDDT-defined disorder, all canonicals with cohort labels"),
        ("Table S7", "Deletion windows and cleavage-site overlap in the emergence cohort"),
        ("Table S8", "STRING protein-protein interactions"),
        ("Table S9", "Gene Ontology and pathway enrichment for the emergence and complete-ablation cohorts against a custom background of all canonical genes. The 36-gene subset of emergence events reproducing on ipSAE was tested under identical settings and returned no significant terms, so it contributes no rows."),
        ("Table S10", "Agreement between the four interface metrics and class transitions"),
        ("Table S11", "Matched controls: random-window null, disorder strata, interface retention, bootstrap intervals"),
        ("Table S12", "Shared splicing events across paralogues, and topology under each assignment"),
        ("Table S13", "Deposited-structure state of the canonical proteins: largest number of copies in any deposited biological assembly, for every canonical carrying a Protein Data Bank entry."),
    ])
    lines = [["Supplementary Information"], [],
             ["Title:", "Structural grammar of truncation across the human homodimer landscape"],
             ["Authors:", "Taner Karagöl, Alper Karagöl"], [], ["Table of Contents"]]
    for k, v in toc.items():
        lines.append([k, v])
    lines += [[], ["Supplementary Figures"],
              ["Figures S1 and S2", "Supplied in the Supplementary Information document, not in this workbook."]]
    lines += [[], ["Notes"],
              ["Table S1", "deletion_window and cleavage_site_overlap are stored as text. A window spanning more than one interval is given as semicolon-separated ranges in canonical numbering."],
              ["Table S1", "topology_alignment is derived by affine-gap global alignment with free end gaps, with a five-residue terminal tolerance."],
              ["Table S1", "canon_ipSAE_min and iso_ipSAE_min are the lower of the two chain directions, matching the release's own definition."],
              ["Table S5", "reproduces_on_ipSAE is TRUE where the canonical falls below 0.6 and the isoform reaches 0.6 or above on ipSAE_min."],
              ["Table S6", "disorder_fraction_pLDDT50 is the fraction of residues with pLDDT below 50, read from the AFDB prediction metadata."]]
    for line in lines:
        ws.append(line)
        for c in ws[ws.max_row]:
            c.font = FONT

    wb.save(DST)
    print("sheets written: %s" % ", ".join(wb.sheetnames))
    for name in wb.sheetnames:
        print("  %-12s %d rows" % (name, wb[name].max_row))


if __name__ == "__main__":
    main()
