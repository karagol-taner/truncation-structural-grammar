"""
12_figure1_landscape.py
=======================
Figure 1 -- ipTM Ablation and Unmasking Landscape.

A scatter of canonical vs. truncated-isoform ipTM, coloured by the
stratified structural outcome.  Points above the y=x line that cross
the 0.6 ipTM threshold are recoloured purple to highlight Gain of
Function (Interface Unmasking) events.

Reads:  AFDB_Truncation_Stratified_Analysis_FIXED.csv
Writes: ipTM_Ablation_and_Unmasking_Landscape.png  (also .svg)

Optional shell prerequisites for headless figure export are listed at
the bottom of the file.

Run order: 12 of 14
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns


def make_landscape(output_basename="ipTM_Ablation_and_Unmasking_Landscape"):
    BASE_PATH = Path('/content/drive/MyDrive/MIT/All Project Files/AFDB dimer')
    if BASE_PATH.exists():
        os.chdir(BASE_PATH)

    print("Loading dataset and identifying Gain of Function (Unmasking)...")
    df = pl.read_csv("AFDB_Truncation_Stratified_Analysis_FIXED.csv")

    # Recolour upward outliers as Gain of Function
    df = df.with_columns(
        pl.when((pl.col("iso_ipTM") > pl.col("canon_ipTM")) &
                (pl.col("iso_ipTM") >= 0.6))
          .then(pl.lit("Gain of Function (Unmasking)"))
          .otherwise(pl.col("structural_outcome"))
          .alias("structural_outcome")
    )
    df_pd = df.to_pandas()

    plt.figure(figsize=(12, 10))
    sns.set_theme(style="whitegrid")

    scientific_palette = {
        "Gain of Function (Unmasking)":         "#6a3d9a",
        "Conserved at >0.8":                    "#08519c",
        "Conserved at >0.7":                    "#4292c6",
        "Conserved at >0.6":                    "#9ecae1",
        "Ablation from >0.8 Canonical":         "#cb181d",
        "Ablation from >0.7 Canonical":         "#fb6a4a",
        "Ablation from >0.6 Canonical":         "#fdd0a2",
        "Monomer Baseline / Stochastic Signal": "#bdbdbd",
    }

    hue_order = [
        "Gain of Function (Unmasking)",
        "Ablation from >0.8 Canonical",
        "Ablation from >0.7 Canonical",
        "Ablation from >0.6 Canonical",
        "Conserved at >0.8",
        "Conserved at >0.7",
        "Conserved at >0.6",
        "Monomer Baseline / Stochastic Signal",
    ]

    sns.scatterplot(
        data=df_pd,
        x='canon_ipTM', y='iso_ipTM',
        hue='structural_outcome',
        palette=scientific_palette,
        hue_order=hue_order,
        alpha=0.7, edgecolor='w', linewidth=0.3, s=60,
    )

    threshold = 0.6
    plt.axhline(threshold, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
    plt.axvline(threshold, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
    plt.plot([0, 1], [0, 1], color='#525252', linestyle=':', alpha=0.6, linewidth=1)

    plt.xlabel('Canonical Homodimer ipTM',          fontsize=12)
    plt.ylabel('Truncated Isoform Homodimer ipTM',  fontsize=12)
    plt.xlim(0, 1)
    plt.ylim(0, 1)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
               frameon=False, fontsize=10)
    plt.tight_layout()

    # Save both PNG (high-res raster) and SVG (vector) outputs
    plt.savefig(f"{output_basename}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_basename}.svg", format='svg', bbox_inches='tight')
    print(f"Saved {output_basename}.png and {output_basename}.svg")
    plt.show()


# ---------------------------------------------------------------------------
# Optional shell prerequisites for headless figure export (Colab):
#   !pip install -U kaleido Plotly
#   !pip install pandas numpy matplotlib seaborn
#   !plotly_get_chrome
#   !sudo apt update && sudo apt-get install \
#         libnss3 libatk-bridge2.0-0 libcups2 libxcomposite1 libxdamage1 \
#         libxfixes3 libxrandr2 libgbm1 libxkbcommon0 libpango-1.0-0 \
#         libcairo2 libasound2
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    make_landscape()
