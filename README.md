# Kunitz HMM Project

This repository contains the workflow, scripts, notebooks, input files, intermediate files, and results for a bioinformatics project focused on the **Kunitz/BPTI protease inhibitor domain**.

The aim of the project was to build a profile Hidden Markov Model (HMM) for the Kunitz domain and evaluate its ability to identify Kunitz-containing proteins in UniProt/Swiss-Prot datasets.

---

## Repository Structure

```text
kunitz-hmm-project/
├── scripts/
│   ├── performance.py
│   └── kunitz_pipeline.sh
│
├── notebooks/
│   ├── data_preparation.ipynb
│   └── roc_plot_final.ipynb
│
├── data/
│   ├── raw/
│   └── processed/
│
├── results/
│
├── figures/
│
└── report/
```

---

## Project Overview

### 1. Structural data collection and filtering

Kunitz-domain structural data were retrieved and filtered using a PDB custom report. The filtered entries were used as the starting point for selecting representative Kunitz-domain sequences and structures.

### 2. Redundancy reduction

The selected sequences were clustered using **MMseqs2** to reduce redundancy and obtain a representative set of Kunitz-domain sequences.

### 3. Structural alignment

Representative Kunitz-domain structures were aligned using **PDBeFold**. The resulting alignment was used as the basis for profile HMM construction.

### 4. Profile HMM construction

A profile Hidden Markov Model was built using **HMMER**:

```bash
hmmbuild kunitz.hmm fasta.seq
```

### 5. UniProt/Swiss-Prot search

The HMM was searched against positive and negative UniProt/Swiss-Prot datasets using `hmmsearch`.

### 6. Validation dataset preparation

HMMER outputs were converted into labelled validation files with the format:

```text
protein_id    evalue    true_label
```

where:

```text
1 = known Kunitz protein
0 = known non-Kunitz protein
```

### 7. Performance evaluation

The model was evaluated across multiple E-value thresholds. Performance was assessed using confusion matrices and classification metrics such as accuracy, MCC, sensitivity, specificity, false positive rate, precision, and F1 score.

### 8. Visualization

The final visualizations include:

- HMM logo generated with **Skylign**
- Semi-log ROC curve generated in Python

---

## Main Files

### Scripts

- `scripts/kunitz_pipeline.sh`  
  Main Bash workflow for HMM construction, HMMER searches, validation-set preparation, and performance evaluation.

- `scripts/performance.py`  
  Python script for calculating threshold-dependent classification metrics.

### Notebooks

- `notebooks/data_preparation.ipynb`  
  Notebook used for data filtering and preparation.

- `notebooks/roc_plot_final.ipynb`  
  Notebook used to generate ROC-related visualizations.

### Data

- `data/raw/`  
  Original input files, including UniProt/Swiss-Prot FASTA files and the PDB custom report.

- `data/processed/`  
  Processed files generated during the workflow, including cluster outputs, alignments, labelled match files, and validation datasets.

### Results

- `results/`  
  HMMER output files, the final profile HMM, and threshold-performance result files.

### Figures

- `figures/`  
  Report-ready figures, including the Skylign HMM logo and ROC curve.

---

## Tools Used

- HMMER
- MMseqs2
- PDBeFold
- Skylign
- Python 3
- NumPy
- pandas
- matplotlib

---

## Notes

The final report contains the full methodology, validation results, threshold selection, figures, and discussion.

## Academic Context

This project was developed as part of the Laboratory of Bioinformatics I - Module 2 taught by professor E. Capriotti at the University of Bologna. 