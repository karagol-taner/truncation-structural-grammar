"""
04_homodimer_mapping.py
=======================
Strict accession-based mapping of human canonical homodimers to their
naturally truncated isoforms, using the AFDB manifest plus the FASTA
sequence lookup.

Run order: 4 of 14
"""

from difflib import SequenceMatcher

import polars as pl


# ---------------------------------------------------------------------------
# Accounting dictionary (transparent statistics for the manuscript)
# ---------------------------------------------------------------------------
stats_counter = {
    "total_rows_processed": 0,
    "missing_sequence": 0,
    "no_canonical_found": 0,
    "failed_truncation_check": 0,
    "successfully_validated_pairs": 0,
}


# ---------------------------------------------------------------------------
# FASTA helpers (small, kept here because the rest of the file consumes them)
# ---------------------------------------------------------------------------
def build_seq_lookup(fasta_path):
    """Index every >header into RAM for instant resolve_seq() lookups."""
    print("Loading FASTA sequences into RAM. This will take a moment...")
    lookup = {}
    with open(fasta_path, 'r') as f:
        header = ""
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    lookup[header] = "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            lookup[header] = "".join(seq)
    print(f"Success! Indexed {len(lookup)} sequences.")
    return lookup


def resolve_seq(acc, seq_map):
    """Try every plausible AFDB FASTA key for a given UniProt accession."""
    keys = [
        f"AFDB:AF-{acc}-F1",
        f"AFDB:AF-{acc.split('-')[0]}-1-F1",
        f"AFDB:AF-{acc}-1-F1",
    ]
    for k in keys:
        if k in seq_map:
            return seq_map[k]
    return None


def is_true_truncation(c_seq, i_seq, threshold=0.95):
    """The shorter sequence must align as a near-perfect substring of the longer."""
    if len(i_seq) >= len(c_seq):
        return False
    if i_seq in c_seq:
        return True
    s = SequenceMatcher(None, c_seq, i_seq)
    match = s.find_longest_match(0, len(c_seq), 0, len(i_seq))
    return (match.size / len(i_seq)) >= threshold


# ---------------------------------------------------------------------------
# Main mapping engine
# ---------------------------------------------------------------------------
def execute_final_homodimer_mapping(df, seq_map):
    print("Executing Strict Accession-Based Mapping (Correcting Gene Ambiguity)...")

    # 1. PRE-PROCESS: Human only, base accession (strip "-N"), keep highest ipTM
    df_prepared = (
        df.filter(pl.col("taxId") == 9606)
        .with_columns(
            pl.col("uniprotAccession").str.replace(r"-\d+$", "").alias("base_acc")
        )
        .sort("ipTM", descending=True)
        .unique(subset=["uniprotAccession"])
    )

    results = []

    # 2. GROUP BY BASE: never compare across families
    for base_tuple, group in df_prepared.group_by("base_acc"):
        base_acc = (str(base_tuple[0]) if isinstance(base_tuple, tuple)
                    else str(base_tuple))

        members = []
        for row in group.iter_rows(named=True):
            stats_counter["total_rows_processed"] += 1
            seq = resolve_seq(row["uniprotAccession"], seq_map)
            if not seq:
                stats_counter["missing_sequence"] += 1
                continue
            members.append({
                "acc":  row["uniprotAccession"],
                "seq":  seq,
                "ipTM": row["ipTM"],
                "gene": row["gene"],
            })

        if len(members) < 2:
            continue

        # 3. Pick the canonical for this family
        canonical_pool = [
            m for m in members
            if m["acc"] == base_acc or m["acc"].endswith("-1")
        ]
        if not canonical_pool:
            stats_counter["no_canonical_found"] += 1
            continue
        canon = max(canonical_pool, key=lambda m: m["ipTM"])

        # 4. Pairwise test every other member against the canonical
        for iso in members:
            if iso["acc"] == canon["acc"]:
                continue
            if is_true_truncation(canon["seq"], iso["seq"]):
                c_val, i_val = canon["ipTM"], iso["ipTM"]

                # Tiered structural outcome label
                if c_val >= 0.8:
                    outcome = ("Ablation from >0.8 Canonical" if i_val < 0.8
                               else "Conserved at >0.8")
                elif c_val >= 0.7:
                    outcome = ("Ablation from >0.7 Canonical" if i_val < 0.7
                               else "Conserved at >0.7")
                elif c_val >= 0.6:
                    outcome = ("Ablation from >0.6 Canonical" if i_val < 0.6
                               else "Conserved at >0.6")
                else:
                    outcome = "Monomer Baseline / Stochastic Signal"

                results.append({
                    "gene":               canon["gene"],
                    "canonical_acc":      canon["acc"],
                    "isoform_acc":        iso["acc"],
                    "canon_ipTM":         c_val,
                    "iso_ipTM":           i_val,
                    "canon_len":          len(canon["seq"]),
                    "iso_len":            len(iso["seq"]),
                    "structural_outcome": outcome,
                })
                stats_counter["successfully_validated_pairs"] += 1
            else:
                stats_counter["failed_truncation_check"] += 1

    return pl.DataFrame(results)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Loading CSV manifest...")
    df_manifest = pl.read_csv("model_entity_metadata_mapping.csv")

    seq_lookup = build_seq_lookup("sequences.fasta")
    final_analysis_df = execute_final_homodimer_mapping(df_manifest, seq_lookup)

    output_filename = "AFDB_Truncation_Stratified_Analysis_FIXED.csv"
    final_analysis_df.write_csv(output_filename)
    print(f"\nSaved clean results to {output_filename}")

    print("\n--- FINAL DATA ACCOUNTING ---")
    for key, val in stats_counter.items():
        print(f"{key.replace('_', ' ').title()}: {val}")
