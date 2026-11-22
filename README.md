# Structural grammar of truncation across the human homodimer landscape

**Taner Karagöl<sup>1,†,\*,‡</sup>, Alper Karagöl<sup>1,¶,\*,‡</sup>**

<sup>1</sup> Istanbul University Istanbul Medical Faculty, Istanbul, Turkey
<sup>†</sup> Present Address: McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA
<sup>¶</sup> Present Address: Broad Institute of MIT and Harvard, Cambridge, MA, USA; Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA
<sup>\*</sup> These authors contribute equally
<sup>‡</sup> To whom the correspondence should be addressed

## Abstract

Alternative splicing and proteolytic truncation generate tens of thousands of protein isoforms in the human proteome. Their consequences for quaternary state, the level at which most protein function operates, have been examined largely one molecule at a time. Using the AlphaFold Protein Structure Database expansion to predicted human homodimers, we compared 5,168 canonical-versus-truncated pairs on four interface metrics distributed with that release. Truncation partitions high-confidence canonical homodimers into structural conservation (56.4% of pairs), complete interface ablation (26.1%) and partial destabilization (17.5%). The partition holds on predicted interface size and not only on interface confidence. Conservation events retain a median 0.938 of the canonical interface residues, destabilization events 0.622, and ablation events 0.000. Interface loss is graded by deletion magnitude, falling from 0.997 of the interface retained below a tenth of the chain removed to 0.192 above a half (Spearman rho = -0.552, p = 4.5 × 10⁻⁶¹). Deletion topology, derived by global alignment, is 53.3% N-terminal and 42.7% C-terminal. C-terminal truncation lowers interface conservation from 65.4% to 37.5%, an asymmetry confined to deletions above the cohort median of 89 residues. Of 208 pairs in which ipTM reports emergence from a sub-threshold canonical baseline, 42 reproduce on ipSAE, the metric on which the release defines its tiers. Canonical disorder content shows no detectable difference between emergence events and controls matched on canonical confidence. Predicted confidence tracks deposited structures, with two or more copies in a deposited assembly for 95.3% of high-confidence canonicals against 49.9% of low-confidence canonicals. Emergence-cohort canonicals are enriched for established homomers, 75.0% against 48.7% for other canonicals of the same confidence, so emergence marks full-length predictions that fail and not interfaces created by truncation. Truncation-associated remodeling of homodimer architecture is a graded, topology-asymmetric property of the human proteome, readable from interface size and operating across gene-functional categories.

**Keywords:** Truncated isoforms, Protein quaternary structure, Intrinsic disorder, Alternative splicing, Homodimer architecture

## Repository layout

```
GitHub/
├── README.md                             This file
├── Supplementary_information.pdf         Supplementary Figures S1 and S2
├── Supplementary_tables.xlsx             Supplementary Tables S1 to S13
├── code/                                 Numbered Python pipeline (run in order) from the first
│                                         release, which named the emergence class "unmasking"
│   ├── 01_setup_environment.py           Working directory and paths
│   ├── 02_data_download.py               AFDB FTP probe, manifest and FASTA download
│   ├── 03_afdb_api_fetcher.py            Per-protein AFDB model fetcher
│   ├── 04_homodimer_mapping.py           Canonical-isoform mapping, writes Full_dataset.csv
│   ├── 05_state_transition_profiling.py  Tier transitions, topology, Tables 1 to 4
│   ├── 06_manuscript_qa_summary.py       Summary and simplified table
│   ├── 07_extract_unmasking.py           Emergence cohort
│   ├── 08_idr_binding_profiling.py       AIUPred and AlphaFold pLDDT disorder profiling
│   ├── 09_domain_annotation.py           UniProt feature mining
│   ├── 10_go_enrichment.py               Functional enrichment
│   ├── 11_go_diagnostic.py               Power sweep and diagnostic enrichment
│   ├── 12_figure1_landscape.py           ipTM landscape panel (auxiliary/)
│   ├── 13_figure2_panels.py              Residues lost, topology, fold enrichment panels
│   └── 14_figure3_panels.py              Topology, disorder and feature-deletion panels
├── revision/                             Multi-metric re-analysis
│   ├── input/                            Cohort tables the re-analysis is checked against
│   │   ├── pairs_master.csv              All 5,168 pairs with outcome labels
│   │   ├── topology_outcome_counts.csv   Topology against outcome
│   │   ├── ablation_cohort.csv           Complete-ablation events
│   │   ├── destabilization_leads.csv     Partial-destabilization leads
│   │   ├── emergence_cohort.csv          Emergence events
│   │   ├── emergence_deletion_windows.csv  Deletion windows and feature overlap
│   │   ├── disorder_scores.csv           AIUPred and pLDDT disorder scores
│   │   └── target_accessions.txt         Accessions to fetch from UniProt and AFDB
│   ├── code/
│   │   ├── fetch_metadata.py             Human subset of the release metadata
│   │   ├── fetch_uniprot.py              Sequences and cleavage-site features
│   │   ├── fetch_plddt.py                pLDDT disorder fractions
│   │   ├── fetch_go.py                   GO annotations
│   │   ├── r01_topology_arbitration.py   Topology assignments against an external criterion
│   │   ├── r02_align_topology.py         Alignment-derived deletion windows
│   │   ├── r03_rule1_recheck.py          Topology against outcome
│   │   ├── r04_crosschecks.py            Outcome labels and window comparison
│   │   ├── r05_ipsae_reclassification.py Classification on ipSAE and interface counts
│   │   ├── r06_landscape_and_lengthmodel.py  Four-class landscape, length-dilution model
│   │   ├── r07_random_window_control.py  Matched null for feature deletion
│   │   ├── r08_disorder_matched_control.py   Disorder within matched strata
│   │   ├── r09_shared_splice_events.py   Paralogous splicing events
│   │   ├── r10_multimetric.py            Agreement across four interface metrics
│   │   ├── r11_rheostat_interface.py     Graded interface loss
│   │   ├── r12_sensitivity.py            Sensitivity and bootstrap intervals
│   │   ├── r14_go_enrichment.py          Hypergeometric enrichment
│   │   ├── r17_disorder_by_topology.py   Disorder deletion within each topology
│   │   ├── r18_emergence_controls.py     Threshold margin, cleavage sites, length confound
│   │   ├── r19_enrichment_summary.py     Summary of the g:Profiler exports
│   │   ├── r20_pdb_crosscheck.py         Predicted confidence against deposited structures
│   │   └── r15_build_supplementary.py    Supplementary workbook
│   └── figcode/
│       ├── figstyle.py                   Shared style and palettes
│       ├── fig1_landscape.py             Figure 1
│       ├── fig2_rheostat.py              Figure 2
│       ├── fig3_topology.py              Figure 3
│       ├── fig4_emergence.py             Figure 4
│       ├── fig5_prepare_models.py        Figure 5 common frame and scale
│       ├── fig5_render_pymol.py          Figure 5 structure panels
│       ├── fig5_render_3dmol.py          Figure 5 structure panels, 3Dmol.js alternative
│       ├── fig5_assemble.py              Figure 5
│       ├── fig6_enrichment.py            Figure 6
│       └── figS1_controls.py             Figure S1
├── dataset/                              Curated input cohorts
│   ├── Full_dataset.csv                  5,168 canonical-versus-truncated pairs
│   ├── Dimer_Gain_Cohort.csv             Emergence cohort
│   └── Ablation_Cohort.csv               Complete-ablation cohort
├── analysis/                             Manuscript tables; Table1 to Table5 and their
│                                         structural_outcome values are from the first release
│   ├── Table1_State_Transition_Matrix.csv
│   ├── Table2_Topology_vs_Biological_Outcome.csv
│   ├── Table3_Top_Destabilization_Leads.csv
│   ├── Table4_Top_Unmasking_Leads.csv
│   ├── Table5_High_Value_Unmasking_Leads.csv
│   ├── Simplified_Manuscript_Table.csv
│   ├── pairs_ipsae.csv                   Every pair on four interface metrics
│   ├── alignment_topology.csv            Alignment-derived windows and topology
│   ├── plddt_disorder.csv                Canonical pLDDT disorder fractions
│   └── pdb_homomer_state.csv             Copies of each canonical in deposited assemblies
├── auxiliary/                            Intermediate and supporting outputs
│   ├── Domain_Annotations_v2*.csv|.txt   UniProt feature-mining outputs
│   ├── IDR_Comprehensive_Scores.csv      Disorder scores
│   ├── IDR_Comprehensive_Profiling.{pdf,png}
│   ├── IDR_Profiling_summary.txt
│   ├── GO_Enrichment_v2_*.csv|.txt
│   ├── GO_Diagnostic_*.csv|.txt
│   ├── gprofiler/                        g:Profiler queries, background and exports
│   ├── enrichment_summary.json           Terms, query sizes and p values per cohort
│   ├── emergence_controls.json           Threshold margin, cleavage sites, length confound
│   ├── disorder_by_topology.json         Disorder deletion within each topology
│   ├── pdb_crosscheck.json               Homomer rate by predicted confidence tier
│   ├── fig5_structures/                  Predicted dimer models shown in Figure 5
│   ├── fig5_oriented/                    Superposed models and their shared scale
│   ├── fig5_panels/                      Rendered panels, one per model
│   ├── fig5_manifest.json                Cases, deletion windows and metrics
│   ├── fig5_plddt.json                   Mean pLDDT per model
│   ├── topology_arbitration.json         Topology under each assignment
│   ├── alignment_topology_summary.json   Counts behind the alignment assignment
│   ├── rule1_recheck.json                Topology against outcome
│   ├── landscape.json                    Four-class landscape on ipSAE
│   ├── ipsae_summary.json                Classification on ipSAE
│   ├── multimetric.json                  Metric agreement and class transitions
│   ├── rheostat_interface.json           Interface retention by deletion magnitude
│   ├── random_window_control.json        Matched null for feature deletion
│   ├── disorder_control.json             Disorder within matched strata
│   ├── sensitivity.json                  Sensitivity and bootstrap intervals
│   ├── shared_splice_events.json         Paralogous splicing events
│   ├── ipTM_Ablation_and_Unmasking_Landscape.png
│   └── string_*                          STRING interaction network
└── figures/                              Manuscript figures
    ├── Figure1_landscape_two_metrics.{png,pdf,svg}
    ├── Figure2_graded_interface_loss.{png,pdf,svg}
    ├── Figure3_topology.{png,pdf,svg}
    ├── Figure4_emergence_discrimination.{png,pdf,svg}
    ├── Figure5_interface_structures.{png,pdf,svg}
    ├── Figure6_enrichment_contrast.{png,pdf,svg}
    ├── FigureS1_controls.{png,pdf,svg}
    └── FigureS2_string_network.png       STRING network, supplementary
```

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `dataset/Full_dataset.csv` | 5,168 | Canonical-versus-truncated homodimer pairs with `canon_ipTM`, `iso_ipTM`, lengths and structural-outcome label. |
| `analysis/pairs_ipsae.csv` | 5,168 | The same pairs with ipSAE, pDockQ, LIS, inter-chain contacts and interface residue counts for both members. |
| `analysis/alignment_topology.csv` | 5,164 | Alignment-derived deletion windows, junction substitutions, topology and cleavage-site overlap. |
| `analysis/plddt_disorder.csv` | 4,115 | Canonical pLDDT disorder fractions from the AFDB prediction metadata. |
| `dataset/Dimer_Gain_Cohort.csv` | 208 | Emergence cohort, pairs crossing an upward ipTM tier boundary. |
| `dataset/Ablation_Cohort.csv` | 61 | Complete-ablation cohort, High to Low transitions. |
| `analysis/pdb_homomer_state.csv` | 4,115 | Largest number of copies of each canonical in any deposited biological assembly, with the count of Protein Data Bank entries behind it. |

Master pairs are derived from the AlphaFold Database homodimer release metadata by accession-based mapping against the AFDB FASTA, with a subsequence criterion at 0.95 identity. The release distributes ipTM, ipSAE, pDockQ, LIS, inter-chain contact counts and interface residue counts for every predicted complex; `revision/code/fetch_metadata.py` extracts the human subset.

## Reproducing the analysis

Scripts use plain relative file names and read their roots from three environment variables.

| Variable | Points at | Used by |
| --- | --- | --- |
| `AFDB_DIMER_DIR` | A working directory for the downloaded AFDB manifest and FASTA | `code/` |
| `AFDB_REVISION_DIR` | A working directory holding `out/` and `data/` for the re-analysis | `revision/code/` |
| `AFDB_INPUT_DIR` | `revision/input` in this repository | `revision/code/` |
| `AFDB_PDB_CACHE` | Where the SIFTS and PDBe downloads are cached, `pdbcache` under `AFDB_REVISION_DIR` by default | `r20_pdb_crosscheck.py` |

`AFDB_REPO_DIR` is optional: the figure scripts and `code/12` to `code/14` resolve the repository root from their own location, and read their inputs from `analysis/`, `dataset/` and `auxiliary/`.

```bash
pip install polars pandas numpy matplotlib seaborn scipy requests \
            biopython openpyxl statsmodels gprofiler-official \
            pillow playwright
playwright install chromium

export AFDB_DIMER_DIR=/path/to/working/directory
export AFDB_REVISION_DIR=/path/to/revision/directory
export AFDB_INPUT_DIR="$PWD/revision/input"

python code/01_setup_environment.py
python code/02_data_download.py
python code/03_afdb_api_fetcher.py
python code/04_homodimer_mapping.py
python code/05_state_transition_profiling.py
python code/06_manuscript_qa_summary.py
python code/07_extract_unmasking.py
python code/08_idr_binding_profiling.py
python code/09_domain_annotation.py
python code/10_go_enrichment.py
python code/11_go_diagnostic.py

python revision/code/fetch_metadata.py
python revision/code/fetch_uniprot.py
python revision/code/fetch_plddt.py
python revision/code/fetch_go.py
python revision/code/r02_align_topology.py
python revision/code/r01_topology_arbitration.py
python revision/code/r03_rule1_recheck.py
python revision/code/r04_crosschecks.py
python revision/code/r05_ipsae_reclassification.py
python revision/code/r06_landscape_and_lengthmodel.py
python revision/code/r07_random_window_control.py
python revision/code/r08_disorder_matched_control.py
python revision/code/r09_shared_splice_events.py
python revision/code/r10_multimetric.py
python revision/code/r11_rheostat_interface.py
python revision/code/r12_sensitivity.py
python revision/code/r14_go_enrichment.py
python revision/code/r17_disorder_by_topology.py
python revision/code/r18_emergence_controls.py
python revision/code/r19_enrichment_summary.py
python revision/code/r20_pdb_crosscheck.py
python revision/code/r15_build_supplementary.py "$AFDB_REVISION_DIR/out" Supplementary_tables.xlsx Supplementary_tables.xlsx

python revision/figcode/fig1_landscape.py
python revision/figcode/fig2_rheostat.py
python revision/figcode/fig3_topology.py
python revision/figcode/fig4_emergence.py
python revision/figcode/fig5_prepare_models.py
python revision/figcode/fig5_render_pymol.py      # or: python revision/figcode/fig5_render_3dmol.py
python revision/figcode/fig5_assemble.py
python revision/figcode/fig6_enrichment.py
python revision/figcode/figS1_controls.py
```

Each script's docstring lists its inputs and outputs. The published Figure 5 panels come from `fig5_render_pymol.py`, which ray-traces each panel in a PyMOL session. `fig5_render_3dmol.py` produces the same panels with 3Dmol.js in headless Chromium and loads the library from a pinned CDN URL; set `MOL3D_JS` to a local copy to run offline.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Data wrangling | `polars`, `pandas`, `numpy` |
| Plotting | `matplotlib`, `seaborn`, `pillow` |
| Statistics | `scipy`, `statsmodels` |
| Structural and sequence | `biopython` (`Bio.Align.PairwiseAligner`, `Bio.PDB.MMCIFParser`) |
| Structure rendering | `pymol`; `playwright` with Chromium and 3Dmol.js for the alternative renderer |
| Spreadsheets | `openpyxl` |
| External APIs | `requests`, `gprofiler-official`, `aria2c` (system), `ftplib` |
| Optional | `AIUPred` for the binding-propensity profiling in `code/08` |

Tested on Python 3.10 and 3.11.

## Data availability

* AlphaFold Database homodimer release metadata and `sequences.fasta` are publicly available from EBI: `ftp://ftp.ebi.ac.uk/pub/databases/alphafold/`.
* UniProtKB feature annotations, GO annotations and FASTA are queried through the UniProt REST API (`https://rest.uniprot.org/`).
* Functional enrichment uses g:Profiler (`https://biit.cs.ut.ee/gprofiler/`) and a hypergeometric test against UniProt GO annotations. The submitted gene lists, the custom background and the exported term tables are under `auxiliary/gprofiler/`.
* Deposited structures are taken from the Protein Data Bank through PDBe (`https://www.ebi.ac.uk/pdbe/`), with chain-level mappings from SIFTS (`https://ftp.ebi.ac.uk/pub/databases/msd/sifts/`).
* The interaction network was generated with STRING v12 (`https://string-db.org/`).
* All cohort files and analysis tables required to reproduce the figures are included under `dataset/`, `analysis/`, `auxiliary/` and `revision/input/`.

## Citation

If you use this code or dataset, please cite:

> Karagöl T., Karagöl A. *Structural grammar of truncation across the human homodimer landscape.* (2026).

## Contact

Questions and issues are welcome. Please open a GitHub Issue or contact the corresponding authors.

## License

The code in this repository is released under the MIT License. Curated datasets and figures are released under CC-BY 4.0. Underlying AFDB, UniProt, g:Profiler and STRING data are subject to the licenses of those upstream providers.
