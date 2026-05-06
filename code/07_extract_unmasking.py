"""
07_extract_unmasking.py
=======================


Run order: 7 of 14
"""

import os
from pathlib import Path

import polars as pl


def extract_unmasking_instances(file_path):
    print(f"Loading {file_path} to isolate Interface Unmasking (Dimer Gain) events...")
    df = pl.read_csv(file_path)

    unmasking_df = df.filter(
        # Low -> Moderate / High
        ((pl.col("canon_ipTM") < 0.6) & (pl.col("iso_ipTM") >= 0.6)) |
        # Moderate -> High
        ((pl.col("canon_ipTM") < 0.8) & (pl.col("iso_ipTM") >= 0.8))
    )

    unmasking_df = unmasking_df.with_columns([
        (pl.col("iso_ipTM") - pl.col("canon_ipTM")).alias("ipTM_gain"),
        (pl.col("canon_len") - pl.col("iso_len")).alias("residues_removed"),
    ])

    unmasking_df = unmasking_df.sort("ipTM_gain", descending=True)
    return unmasking_df


if __name__ == "__main__":
    BASE_PATH = Path('/content/drive/MyDrive/MIT/All Project Files/AFDB dimer')
    os.chdir(BASE_PATH)

    unmasking_results = extract_unmasking_instances(
        "AFDB_Truncation_Stratified_Analysis_FIXED.csv"
    )

    unmasking_results.write_csv("Table5_High_Value_Unmasking_Leads.csv")

    print(f"\nSuccessfully isolated {len(unmasking_results)} validated instances.")
    print("Top 10 Structural Activation Leads (Interface Unmasking):")

    print(unmasking_results.select([
        "gene", "canonical_acc", "isoform_acc",
        "canon_ipTM", "iso_ipTM", "ipTM_gain", "residues_removed",
    ]).head(10))
