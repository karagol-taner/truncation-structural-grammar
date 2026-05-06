# A structural grammar of truncation across the human homodimer landscape

**Taner Karagöl & Alper Karagöl**


## Abstract

Alternative splicing and proteolytic truncation generate tens of thousands of protein isoforms in the human proteome, but the structural consequences for quaternary state, the level at which most signaling, enzymatic and regulatory function operates, have largely been examined one molecule at a time. Leveraging the recent expansion of the AlphaFold Database to predicted human homodimers, we systematically compared 5,168 canonical-versus-truncated homodimer pairs across the human proteome. In high-confidence canonical homodimers, truncation is associated with predicted structural conservation in 56.4% of pairs (mean 85 residues lost), complete interface ablation in 26.1% (mean 178 residues lost), and partial destabilization in 17.5% (mean 134 residues lost); a distinct fourth class (4.0% of the dataset, n = 208) shows truncation-associated emergence of a predicted high-confidence interface from a sub-threshold canonical baseline. Two reproducible rules govern these transitions: a topological asymmetry in which N-terminal losses are preferentially enriched ~1.6-fold in interface preservation while C-terminal losses are rare overall (~6% of pairs) and modestly under-represented in the conservation class, and a biophysical rule in which emergence-class proteins show substantially elevated intrinsic disorder content relative to ablation-class proteins, as measured by both AlphaFold pLDDT-defined disorder of the canonical structure (Cohen's d ≈ 1.39) and AIUPred peak binding propensity of the truncated isoform (Cohen's d ≈ 0.65). Formal pathway enrichment recovered only a small nucleotide-metabolism signal, indicating that these rules operate across diverse gene-functional categories. Truncation-associated remodeling of homodimer architecture thus constitutes a structural grammar of the human proteome rather than a specialty of any single regulatory family.

---

## Repository layout

```
GitHub/
├── README.md                        This file
├── Supplementary_Tables.xlsx        Supplementary tables for the manuscript
├── code/                            Numbered Python pipeline (run in order)
│   ├── 01_setup_environment.py      Environment / Drive / paths
│   ├── 02_data_download.py          AFDB FTP probe + manifest + FASTA download
│   ├── 03_afdb_api_fetcher.py       Per-protein AFDB model fetcher
│   ├── 04_homodimer_mapping.py      Strict canonical-isoform mapping
│   ├── 05_state_transition_profiling.py   Tier transitions, topology, Tables 1-4
│   ├── 06_manuscript_qa_summary.py  Q&A summary + simplified table
│   ├── 07_extract_unmasking.py      Dimer-Gain (Unmasking) cohort -> Table 5
│   ├── 08_idr_binding_profiling.py  AIUPred + AlphaFold pLDDT pipeline
│   ├── 09_domain_annotation.py      UniProt feature mining + alignment-based deletion windows
│   ├── 10_go_enrichment.py          g:Profiler primary + low-tier robustness
│   ├── 11_go_diagnostic.py          Power sweep + diagnostic enrichment
│   ├── 12_figure1_landscape.py      Figure 1 ipTM landscape
│   ├── 13_figure2_panels.py         Figure 2 (residues / topology / forest)
│   ├── 14_figure3_panels.py         Figure 3 (5-panel emergence cohort)
│   └── AFDB_truncated_dimer.ipynb   Original Colab notebook
├── dataset/                         Curated input cohorts
│   ├── Full_dataset.csv             5,168 canonical-vs-truncated pairs (master)
│   ├── Dimer_Gain_Cohort.csv        Emergence (Interface Unmasking) cohort
│   └── Ablation_Cohort.csv          Complete Interface Ablation cohort
├── analysis/                        Manuscript tables (results)
│   ├── Table1_State_Transition_Matrix.csv
│   ├── Table2_Topology_vs_Biological_Outcome.csv
│   ├── Table3_Top_Destabilization_Leads.csv
│   ├── Table4_Top_Unmasking_Leads.csv
│   ├── Table5_High_Value_Unmasking_Leads.csv
│   └── Simplified_Manuscript_Table.csv
├── auxiliary/                       Intermediate and supporting outputs
│   ├── Domain_Annotations_v2*.csv|.txt    UniProt feature-mining outputs
│   ├── IDR_Comprehensive_Scores.csv       AIUPred + AlphaFold disorder scores
│   ├── IDR_Comprehensive_Profiling.{pdf,png}
│   ├── IDR_Profiling_summary.txt          Welch / Mann-Whitney / Cohen's d
│   ├── GO_Enrichment_v2_*.csv|.txt        Primary + low-tier robust enrichment
│   ├── GO_Diagnostic_*.csv|.txt           Power sweep + diagnostic enrichment
│   ├── ipTM_Ablation_and_Unmasking_Landscape.png
│   └── string_*                            STRING physical-interaction network
└── figures/                         Final manuscript figures (PNG)
    ├── Figure1_big.png              ipTM landscape scatter
    ├── Figure2_big.png              Residues / topology / forest panels
    ├── Figure3_big.png              Emergence-cohort biophysical panels
    └── Figure4.png                  STRING physical-interaction network
```

---

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `dataset/Full_dataset.csv` | 5,168 | Validated canonical–truncated homodimer pairs with `canon_ipTM`, `iso_ipTM`, lengths, structural-outcome label. |
| `dataset/Dimer_Gain_Cohort.csv` | 208 | Emergence (Interface-Unmasking) cohort — subset where the truncated isoform crosses an upward ipTM tier boundary. |
| `dataset/Ablation_Cohort.csv` | – | Complete Interface Ablation cohort — High → Low transitions. |

Master pairs are derived from the AlphaFold Database 2026 multimer manifest (`model_entity_metadata_mapping.csv`, ~4.6 GB, ~31 M rows) by strict accession-based mapping against the full AFDB FASTA (`sequences.fasta`, ~118 GB) and SequenceMatcher truncation validation (≥ 0.95 substring identity). Detailed accounting is printed by `code/04_homodimer_mapping.py`.

---

## Reproducing the analysis

All scripts use plain relative file names and assume a single working directory.

The notebook `code/AFDB_truncated_dimer.ipynb` contains the same cells in a Colab-ready form.

---

## Software dependencies

| Domain | Packages |
| --- | --- |
| Data wrangling | `polars`, `pandas`, `numpy` |
| Plotting | `matplotlib`, `seaborn` |
| Statistics | `scipy`, `statsmodels` |
| Structural / sequence | `biopython` (pairwise2), `metapredict`, `AIUPred` |
| External APIs | `requests`, `aria2c` (system), `ftplib` |
| Functional enrichment | `gprofiler-official` |

Tested on Python 3.10 (Colab) and Python 3.11.

---

## Data availability

* AlphaFold Database multimer manifest and `sequences.fasta` are publicly available from EBI: `ftp://ftp.ebi.ac.uk/pub/databases/alphafold/`.
* UniProtKB feature annotations and FASTA are queried via the UniProt REST API (`https://rest.uniprot.org/`).
* Functional enrichment uses g:Profiler (`https://biit.cs.ut.ee/gprofiler/`).
* Physical-interaction network was generated with STRING v12 (`https://string-db.org/`).
* All curated cohort files and analysis tables required to reproduce the figures are included under `dataset/` and `analysis/`.

---

## Citation

Karagöl T., Karagöl A. *A structural grammar of truncation across the human homodimer landscape.* (2026).


---

## License

The code in this repository is released under the MIT License. Underlying AFDB, UniProt, g:Profiler, and STRING data are subject to the licenses of those upstream providers.
