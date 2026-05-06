"""
05_state_transition_profiling.py
================================
Granular state-transition profiling of the validated truncation pairs
produced in 04_homodimer_mapping.py.

For every (canonical, isoform) pair this step computes:
  * delta_ipTM, ipTM_gain and delta_len
  * canon_tier / iso_tier (High >= 0.8, Moderate 0.6-0.79, Low < 0.6)
  * tier_transition  (e.g. "High -> Moderate")
  * biological_outcome (Conservation, Destabilization, Ablation,
    Gain-of-Function / Interface Unmasking, Baseline)
  * truncation_topology (N-terminal loss / Internal / C-terminal loss)
    derived directly from the FASTA sequences.

The script writes the four manuscript tables consumed by the figure code:
  Table1_State_Transition_Matrix.csv
  Table2_Topology_vs_Biological_Outcome.csv
  Table3_Top_Destabilization_Leads.csv
  Table4_Top_Unmasking_Leads.csv

Run order: 5 of 14
"""

import polars as pl


# ---------------------------------------------------------------------------
# Small FASTA helpers (kept local so the file is self-contained)
# ---------------------------------------------------------------------------
def build_seq_lookup(fasta_path):
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
    keys = [
        f"AFDB:AF-{acc}-F1",
        f"AFDB:AF-{acc.split('-')[0]}-1-F1",
        f"AFDB:AF-{acc}-1-F1",
    ]
    for k in keys:
        if k in seq_map:
            return seq_map[k]
    return None


# ---------------------------------------------------------------------------
# Main profiling function
# ---------------------------------------------------------------------------
def perform_state_transition_profiling(df, seq_map):
    print("Initiating Granular State Transition Profiling...")

    # 1. Delta metrics & dimer-gain
    df = df.with_columns([
        (pl.col("canon_ipTM") - pl.col("iso_ipTM")).alias("delta_ipTM"),
        (pl.col("iso_ipTM")   - pl.col("canon_ipTM")).alias("ipTM_gain"),
        (pl.col("canon_len")  - pl.col("iso_len")).alias("delta_len"),
    ])

    # 2. ipTM tier classification
    df = df.with_columns([
        pl.when(pl.col("canon_ipTM") >= 0.8).then(pl.lit("High"))
          .when(pl.col("canon_ipTM") >= 0.6).then(pl.lit("Moderate"))
          .otherwise(pl.lit("Low")).alias("canon_tier"),

        pl.when(pl.col("iso_ipTM") >= 0.8).then(pl.lit("High"))
          .when(pl.col("iso_ipTM") >= 0.6).then(pl.lit("Moderate"))
          .otherwise(pl.lit("Low")).alias("iso_tier"),
    ])

    # 3. Transition labels
    df = df.with_columns(
        pl.concat_str([pl.col("canon_tier"), pl.lit(" -> "), pl.col("iso_tier")])
          .alias("tier_transition")
    )

    # 4. Biological outcome
    df = df.with_columns(
        pl.when(pl.col("tier_transition") == "High -> High").then(pl.lit("Conserved High-Confidence Dimer"))
          .when(pl.col("tier_transition") == "High -> Moderate").then(pl.lit("Interface Destabilization (Attenuated)"))
          .when(pl.col("tier_transition") == "High -> Low").then(pl.lit("Complete Interface Ablation"))
          .when(pl.col("tier_transition") == "Moderate -> Moderate").then(pl.lit("Conserved Moderate Dimer"))
          .when(pl.col("tier_transition") == "Moderate -> Low").then(pl.lit("Ablation from Moderate"))
          .when(pl.col("tier_transition").is_in(["Low -> High", "Low -> Moderate", "Moderate -> High"]))
            .then(pl.lit("Gain of Function (Interface Unmasking)"))
          .otherwise(pl.lit("Non-Dimer Baseline (Low -> Low)"))
          .alias("biological_outcome")
    )

    # 5. Topology (where is the cut?)
    def get_topology(struct):
        c_seq = resolve_seq(struct["canonical_acc"], seq_map)
        i_seq = resolve_seq(struct["isoform_acc"],  seq_map)
        if not c_seq or not i_seq:
            return "Unknown Topology"
        if c_seq.startswith(i_seq):
            return "C-Terminal Loss"
        elif c_seq.endswith(i_seq):
            return "N-Terminal Loss"
        else:
            return "Internal Deletion / Splice"

    df = df.with_columns(
        pl.struct(["canonical_acc", "isoform_acc"])
          .map_elements(get_topology, return_dtype=pl.String)
          .alias("truncation_topology")
    )

    # 6. Manuscript Table 1 -- State transition matrix
    transition_matrix = df.group_by(
        ["canon_tier", "tier_transition", "biological_outcome"]
    ).agg([
        pl.len().alias("pair_count"),
        pl.col("delta_ipTM").mean().round(3).alias("avg_ipTM_drop"),
        pl.col("ipTM_gain").mean().round(3).alias("avg_ipTM_gain"),
        pl.col("delta_len").mean().round(1).alias("avg_residues_lost"),
    ]).sort(["canon_tier", "pair_count"], descending=[False, True])

    tier_totals = df.group_by("canon_tier").agg(pl.len().alias("starting_tier_total"))
    transition_matrix = transition_matrix.join(tier_totals, on="canon_tier")
    transition_matrix = transition_matrix.with_columns(
        ((pl.col("pair_count") / pl.col("starting_tier_total")) * 100)
        .round(2).alias("percent_of_starting_tier")
    ).select([
        "canon_tier", "tier_transition", "biological_outcome",
        "pair_count", "starting_tier_total", "percent_of_starting_tier",
        "avg_ipTM_drop", "avg_ipTM_gain", "avg_residues_lost",
    ])

    # 7. Manuscript Table 2 -- Topology vs biological outcome
    topo_analysis = df.group_by(
        ["truncation_topology", "biological_outcome"]
    ).agg([pl.len().alias("count")])\
     .sort(["truncation_topology", "count"], descending=[False, True])

    topo_totals = df.group_by("truncation_topology").agg(pl.len().alias("topo_total"))
    topo_analysis = topo_analysis.join(topo_totals, on="truncation_topology")
    topo_analysis = topo_analysis.with_columns(
        ((pl.col("count") / pl.col("topo_total")) * 100)
        .round(2).alias("percent_of_topology")
    ).select([
        "truncation_topology", "biological_outcome",
        "count", "topo_total", "percent_of_topology",
    ])

    # 8. Lead lists for the manuscript
    destabilized_leads = df.filter(
        pl.col("biological_outcome") == "Interface Destabilization (Attenuated)"
    ).sort("canon_ipTM", descending=True)

    unmasking_leads = df.filter(
        pl.col("biological_outcome") == "Gain of Function (Interface Unmasking)"
    ).sort("iso_ipTM", descending=True)

    return df, transition_matrix, topo_analysis, destabilized_leads, unmasking_leads


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    seq_lookup = build_seq_lookup("sequences.fasta")

    print("Loading data for Granular State Transition Analysis...")
    df_stratified = pl.read_csv("AFDB_Truncation_Stratified_Analysis_FIXED.csv")

    full_df, transition_table, topo_table, destab_leads, unmask_leads = (
        perform_state_transition_profiling(df_stratified, seq_lookup)
    )

    transition_table.write_csv("Table1_State_Transition_Matrix.csv")
    topo_table.write_csv("Table2_Topology_vs_Biological_Outcome.csv")
    destab_leads.write_csv("Table3_Top_Destabilization_Leads.csv")
    unmask_leads.write_csv("Table4_Top_Unmasking_Leads.csv")

    print("\n--- TABLE 1: STATE TRANSITION MATRIX (INCLUDING GAIN) ---")
    print(transition_table)

    print("\n--- TABLE 2: HOW TOPOLOGY DICTATES THE OUTCOME ---")
    print(topo_table)
