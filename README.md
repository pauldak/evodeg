# evodeg - QCAP Degron Structural Conservation Pipeline

**Associated publication:**
Shir A. & Ravid T., *"Structural conservation of Degradation Signals is a key determinant of Quality Control-Associated Proteolysis"*, Nature Communications (in review)

A MATLAB pipeline for cross-species structural, sequence, and functional comparison of predicted QCAP degron regions across eukaryotic orthologous protein groups. Using AlphaFold structural models and QCDPred degron predictions, the pipeline quantifies fragment-level conservation across seven species spanning fungal and metazoan lineages, and tests whether degron regions are more conserved than the surrounding protein background.

---

## Quick Start - Run the Demo

No PyMOL or QCDPred required for any demo option.

### Option 1 - Run in browser (no installation required)

Click **Reproducible Run** in the Code Ocean capsule - no account or installation needed.

> **[Open Code Ocean capsule →](https://codeocean.com/capsule/evodeg)**
> *(link will be active upon publication)*

All 6 figures are saved automatically to the results directory.

---

### Option 2 - Run locally in MATLAB

1. Go to the **[Code Ocean capsule](https://codeocean.com/capsule/evodeg)**
2. Click **Download** to get the full capsule as a zip file
3. Extract and open in MATLAB

```matlab
% Figures-only demo - ~10 seconds, Statistics toolbox only
demo

% Full pipeline demo on 2 pre-computed HOGs - ~2–5 min
addpath(genpath(pwd))
demo_pipeline
```

---

### What the Demo Produces

**`demo.m` - 6 figures from a representative pre-computed subset (950 rows per organism):**

| Figure | Content | Statistical test |
|---|---|---|
| 1 | RMSD distributions - degron vs global intervals, per organism (box plots) | Wilcoxon rank-sum |
| 2 | Mean RMSD per organism - paired scatter with SEM | Wilcoxon signed-rank |
| 3 | Sequence alignment score distributions, per organism (box plots) | Wilcoxon rank-sum |
| 4 | Mean alignment score per organism - paired scatter | Wilcoxon signed-rank |
| 5 | Successful vs failed structural coupling - stacked bar chart | Chi-square |
| 6 | Evolutionary rate of change: structural / sequence / functional | Paired t-test + Cohen's d |

**`demo_pipeline.m` - 4 figures freshly computed from 2 HOGs:**
Runs Phases 1–3 of the pipeline live (interval detection, structural coupling, RMSD and alignment scoring) on two pre-computed protein groups, then produces RMSD and alignment score comparison figures. The PyMOL structural alignment step is the only step not reproduced - its outputs are provided pre-computed in `demo_data/hog_data/`. Expected runtime: ~10–15 minutes on a standard desktop (8 cores).

---

## System Requirements

### Software Dependencies

| Dependency | Version tested | Purpose |
|---|---|---|
| MATLAB | R2024a | Core runtime |
| Bioinformatics Toolbox | Required | `fastaread`, `pdbread`, `localalign`, `nwalign` |
| Parallel Computing Toolbox | Required | `parfor` loops in `pdb_superposition` |
| Statistics and Machine Learning Toolbox | Required | `ranksum`, `signrank`, `anova1`, `ttest` |
| Internet access | Required for data acquisition | AlphaFold EBI and UniProt downloads via `blast_cif_downloader.m` |
| PyMOL | 2.x (Schrödinger) | Global and interval-level structural alignment (full pipeline only) |
| QCDPred | - | Degron probability prediction via `prob_calculator_with_fasta_input` (full pipeline only) |

**Tested on:** Windows 10/11, MATLAB R2024a

**Non-standard hardware:** None required. A multi-core CPU is recommended; `parfor` loops utilise all available cores automatically. Full pipeline runtime scales with core count.

**Typical install time:** No installation beyond MATLAB itself is required. MATLAB installation takes approximately 30–60 minutes depending on selected toolboxes.

---

## Repository Structure

```
evodeg/
├── *.m                             Pipeline and analysis functions
├── run                             Code Ocean capsule entry point
├── data/                           Pre-computed tables and reference files
│   ├── hog_protein_id_table.xlsx               Gene group → UniProt ID mapping (228 HOGs × 7 species)
│   ├── evolution_table_column_names.txt     Column names for the evolution comparison table
│   ├── global_coupling_table_degron_alignment.xlsx
│   ├── global_deg_coupling_table_degron_alignment.xlsx
│   ├── global_coupling_table_regular_alignment.xlsx
│   ├── global_deg_coupling_table_regular_alignment.xlsx
│   ├── evolutionary_rates_summary.xlsx
│   ├── pairwise_global_rmsd_matrix.xlsx               7×7 global RMSD matrix (all organism pairs)
│   ├── secondary_structure_degron_matched.xlsx        Secondary structure - matched degron intervals
│   ├── secondary_structure_random_control.xlsx        Secondary structure - randomised control intervals
│   ├── secondary_structure_random_per_organism.xlsx             Per-organism secondary structure counts (random)
│   ├── secondary_structure_degron_per_organism.xlsx             Per-organism secondary structure counts (degron)
│   ├── secondary_structure_random_summary.xlsx              Aggregated secondary structure (random, averaged)
│   └── secondary_structure_degron_summary.xlsx              Aggregated secondary structure (degron, averaged)
├── docs/
│   └── penalty_matrix_for_matlab.xlsx   QCDPred substitution matrix (main pipeline)
├── demo_data/
│   ├── *.xlsx                      Pre-computed coupling tables (sampled subset)
│   ├── penalty_matrix_for_matlab.xlsx
│   └── hog_data/
│       ├── 1377at2759/             HOG 1: pre-computed structural files, 7 organisms
│       └── 604at2759/              HOG 2: pre-computed structural files, 7 organisms
└── multi_evolution_tables/         All-vs-all coupling tables (Figure 5 inputs)
                                    28 files: 7 zero organisms × 2 alignment types ×
                                    2 table types (global / degron). Filename format:
                                    $<zero_organism>$_@<alignment_type>@_<table>.xlsx
```

---

## Setup (Full Pipeline)

1. Clone or download this repository into the MATLAB path
2. Open `addresses.m` and set `main_folder` to the absolute path of the project root on your local machine - all other paths are derived from this value
3. Set `is_exporting_local_files = 1` to write output files to disk, or `= 0` for a dry run
4. Place proteome FASTA files (one per species, downloaded from UniProt) in `proteomes/`
5. Run `addresses()` once to generate `folder_paths.mat` - all pipeline functions load this file at runtime

---

## Required Source Files

The following files must be present before running the full pipeline:

| File | Location | Description |
|---|---|---|
| `data/hog_protein_id_table.xlsx` | `data/` | Maps each HOG ID to one UniProt protein ID per organism (228 rows × 7 species) |
| `docs/penalty_matrix_for_matlab.xlsx` | `docs/` | QCDPred-derived substitution matrix (amino acid order: ARNDCQEGHILKMFPSTWYV) |
| `data/evolution_table_column_names.txt` | `data/` | Column names for the evolution comparison table, one per line |
| `data/global_coupling_table_degron_alignment.xlsx` | `data/` | Pre-computed coupling results - sliding windows, QCDPred alignment |
| `data/global_deg_coupling_table_degron_alignment.xlsx` | `data/` | Pre-computed coupling results - degron windows, QCDPred alignment |

Species proteome FASTA files must be placed in `proteomes/` before running `blast_cif_downloader.m`. AlphaFold PDB and CIF files are downloaded automatically.

---

## Execution Order

The pipeline has four phases. Each depends on the outputs of the previous one.

### Phase 0 - One-time Data Acquisition

```
ortho_db_search_result_extracter   → Download ortholog FASTA files from OrthoDB for 7 species
orthodb_windows_folder_generator   → Create directory structure per HOG and organism
blast_cif_downloader               → For each ortholog: match to species proteome, download
                                     AlphaFold PDB + CIF, UniProt FASTA, run QCDPred
pdbtx_all_protein_initation        → Generate PyMOL script for global alignment of all proteins
                                     to S. cerevisiae; produces _tx.pdb and _tx.aln files
                                     (~2.5 hours for the full 228 HOG dataset)
```

### Phase 1 - Interval Definition and Coupling

```matlab
pdb_superposition()    % Runs first_run() internally; approximately 3 minutes
```

For each HOG (parallelised with `parfor`):
1. `initiate_variables` - Initialise the data structure with global parameters
2. `get_all_crnt_variables` - Load per-organism file paths and protein sequences
3. `get_dist_matrix` - Compute pairwise Cα distance matrix using globally-aligned structures (`_tx.pdb`)
4. `get_global_rmsd_scores` - Extract conserved positions from CLUSTAL alignment (`_tx.aln`); compute whole-protein RMSD
5. `amino_acid_sliding_window` - Generate 32 aa sliding windows (1 aa step) on *S. cerevisiae*
6. `get_degrons_intervals` - Identify degron regions (QCDPred score > 0.85, ±8 aa expansion, contiguous merging)
7. `get_coupling_data` - Map each reference interval to the structurally equivalent region in each comparison species using the Cα distance matrix
8. `create_interval_pymol_files` - Write PyMOL scripts for interval-level alignment; produces `_int_tx.pdb` and `_int_tx.aln`

After parallel loop: `combine_all_pml` → `run_pymol_comnds` (~2.5 hours)

### Phase 2 - Interval-Level RMSD and Sequence Alignment

```matlab
pdb_superposition()    % Runs second_run(hogs_array) internally; approximately 12 minutes
```

For each HOG (parallelised):
1. `dist_matrix_interval_iterator` - Compute Cα distance matrices for each coupled interval pair
2. `get_interval_rmsd_scores('all')` / `('deg')` - Compute interval RMSD and global-projection RMSD
3. `set_alignment_score('all')` / `('deg')` - Compute local sequence alignment scores using QCDPred substitution matrix
4. `export_coupling_results('all')` / `('deg')` - Write per-organism CSV files

### Phase 3 - Table Assembly

```matlab
combine_coupling_table()    % Aggregates all CSVs into global Excel tables
```

Produces:
- `global_coupling_table.xlsx` - sliding-window intervals
- `global_deg_coupling_table.xlsx` - degron intervals

### Phase 4 - Statistical Analysis and Visualization

Run independently in any order after Phase 3:

| Script | Input | Figures produced | Method |
|---|---|---|---|
| `plot_rmsd_alignment_comparisons.m` | Coupling tables | RMSD and alignment box plots + paired scatter (Figures 3A–3D) | Wilcoxon rank-sum + signed-rank |
| `degron_scnd_struct_comparison.m` | Degron coupling table + CIF files | Degron overlap vs random; secondary structure composition | Binomial CI; Chi-square |
| `error_plot_bar_graph.m` | Coupling tables | Stacked bar: match success rate | Chi-square |
| `plot_fig_pie_chart.m` | `evolutionary_rates_summary.xlsx` | Evolutionary slope comparison (Figures 5C–5F) | Paired t-test, Cohen's d, ANOVA |
| `evo_struct_align.m` | `multi_evolution_tables/` | All-vs-all 7×7 distance matrices; divergence trend plots (Figure 5 inputs) | Mean RMSD / alignment score per organism pair |
| `plot_structural_evotime_axis.m` | `pairwise_global_rmsd_matrix.xlsx` | Structure-derived evolutionary time axis (Figure 5B) | Mean RMSD ordering |
| `generate_evo_time_vs_parameters.m` | `hogs_array.mat` | Ribbon plots of conservation vs evolutionary time | Linear regression (slope per HOG) |
| `stats_extract.m` | `evolutionary_rates_summary.xlsx` | Console statistics for Figures 5C–5F | Paired t-test, Cohen's d, ANOVA + Tukey |

---

## Key Parameters

All parameters are set in `initiate_variables.m` and `addresses.m`.

### Pipeline mode (`addresses.m`)

| Parameter | Default | Description |
|---|---|---|
| `is_multi_evo` | 0 | 0 = standard pipeline (*S. cerevisiae* reference); 1 = multi-evolution mode (any reference species) |
| `zero_organism` | `'Saccharomyces_cerevisiae'` | Reference organism. Change for each multi-evolution run, then re-run `addresses()` |
| `is_exporting_local_files` | 0 | 0 = dry run (no files written); 1 = write PDB, ALN, CSV, and PML files to disk |
| `align_to_use` | 2 | Sequence alignment scoring: 1 = standard BLOSUM; 2 = QCDPred-derived custom matrix |
| `full_independant_cycle_run` | 0 | 1 = run `pdbtx_all_protein_initation` + `delete_empty_hogs` at pipeline startup (multi-evolution only) |
| `warning_status` | `'off'` | Warning suppression level during parallel execution |

### Analysis parameters (`initiate_variables.m`)

| Parameter | Value | Description |
|---|---|---|
| `prediction_threshold` | 0.85 | QCDPred logit-smooth score cutoff for degron classification |
| `segment_size` | 32 aa | Sliding window length (based on mean degron size of ~31.5 aa) |
| `gap_coupling` | 8 aa | Expansion applied to each side of a structurally coupled interval boundary |
| `minimal_angstrom_dist` | 2 Å | Minimum Cα distance threshold during coupling search |
| `overlap_value` | 0.25 | Minimum fractional overlap to count a structural interval as matching a degron |

---

## Multi-Evolution Mode

The pipeline supports running with any of the 7 species as the reference (zero) organism, which produces all-vs-all pairwise structural distance matrices used in Figure 5.

**To run one multi-evolution cycle:**

1. In `addresses.m`, set:
   ```matlab
   is_multi_evo = 1;
   zero_organism = 'Homo_sapiens';   % or any other species
   full_independant_cycle_run = 1;   % if PyMOL alignment files need to be generated
   is_exporting_local_files = 1;
   ```
2. Run `addresses()` to update `folder_paths.mat`
3. Run `pdb_superposition()`
4. Repeat for each of the 7 species as reference

Pre-computed coupling tables for all 7 reference organisms (2 alignment types × 2 table types = 28 files total) are provided in `multi_evolution_tables/`. These are the direct inputs to `evo_struct_align.m` for the Figure 5 all-vs-all analysis.

**Organism order note:** When `is_multi_evo = 1`, the pipeline uses `organism_list_evo_ordered` (defined in `addresses.m`), in which *Nematostella vectensis* and *C. elegans* are swapped relative to `organism_list` to reflect structural divergence order from *S. cerevisiae*.

---

## Output Files

| File | Generated by | Description |
|---|---|---|
| `folder_paths.mat` | `addresses.m` | Serialised workspace: all paths and parameters |
| `<id>_tx.pdb` / `_tx.aln` | PyMOL | Globally-aligned structure and CLUSTAL alignment per protein |
| `<id>_<interval>_int_tx.pdb` / `_int_tx.aln` | PyMOL | Interval-aligned structure and alignment |
| `coupling_<id>_all.csv` | `export_coupling_results.m` | Per-organism coupling data (sliding windows) |
| `coupling_<id>_deg.csv` | `export_coupling_results.m` | Per-organism coupling data (degron windows) |
| `global_coupling_table.xlsx` | `combine_coupling_table.m` | Aggregated coupling data - all intervals |
| `global_deg_coupling_table.xlsx` | `combine_coupling_table.m` | Aggregated coupling data - degron intervals |
| `hogs_array.mat` | `pdb_superposition.m` | Full analysis data structure (~70 GB per reference organism) |
| `secondary_structure_degron_matched.xlsx` | `degron_scnd_struct_comparison.m` | Secondary structure of matched degron intervals |
| `secondary_structure_random_control.xlsx` | `degron_scnd_struct_comparison.m` | Secondary structure of length-matched random intervals |
| `statistic_test.txt` | `plot_rmsd_alignment_comparisons.m` | Per-organism Wilcoxon p-values |

---

## Species

| Species | Role |
|---|---|
| *Saccharomyces cerevisiae* | Reference (zero organism) |
| *Amphimedon queenslandica* | Comparison |
| *Nematostella vectensis* | Comparison |
| *Caenorhabditis elegans* | Comparison |
| *Drosophila melanogaster* | Comparison |
| *Danio rerio* | Comparison |
| *Homo sapiens* | Comparison |

**Full dataset:** 228 orthologous protein groups (HOGs), each with one representative protein per species. Full pipeline runtime approximately 3 hours on a standard desktop with 8 cores.

**Reproducibility:** The pipeline is fully deterministic. No fixed random seed is required - the only stochastic component is the randomised control interval generation in `degron_scnd_struct_comparison.m`, which does not affect the primary structural or sequence conservation analyses. Replication of all statistical conclusions is robust across independent runs.

---

## License

MIT License. Copyright 2026 Shir Armony. See `LICENSE` for full terms. The MATLAB runtime is subject to MathWorks licence terms separately.
