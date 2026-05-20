# lsgkm
Sequence-based modelling of H3K27me3 chromatin repression in vascular smooth muscle cells enables prioritisation of cardiovascular disease-associated GWAS variants
Here is a **clean, publication-ready GitHub README.md** tailored to your exact pipeline (LS-GKM + variant prioritisation + TF motif disruption + figures).

You can **copy-paste this directly** into your repo.

***

# H3K27me3 LS‑GKM Variant Prioritisation Pipeline (VSMC)

This repository contains a complete, reproducible pipeline for:

* Training an **LS‑GKM (gkm-SVM)** model on H3K27me3 ChIP‑seq data
* Predicting regulatory impact of **cardiovascular GWAS variants**
* Computing **delta-scores (alt − ref)**
* Performing **TF motif disruption analysis (FIMO + JASPAR)**
* Generating publication-ready **figures (ROC, PR, distributions, prioritisation)**

***

# Overview

We trained a sequence-based classifier to model **H3K27me3 activity in vascular smooth muscle cells (VSMCs)** and used it to:

1. Learn sequence patterns underlying broad H3K27me3 peaks
2. Score **200 bp genomic sequences**
3. Predict the **functional impact of genetic variants**
4. Identify **transcription factor (TF) motif disruptions**

***

# 🔬 Pipeline Summary

## Step 1 — Peak Processing

* Input: nf-core ChIP-seq output
* Consensus peaks filtered using:

```
num_samples ≥ 2 (≥2 of 3 replicates)
```

* Reason: ENCODE standard for broad histone marks
* Result:
  * Total peaks: **8563**
  * High-confidence peaks: **376**
  * Final after blacklist filtering: **375**

***

## Step 2 — Training Data

| Dataset            | Count |
| ------------------ | ----- |
| Positive sequences | 375   |
| Negative sequences | 375   |

* Sequence length: **200 bp**
* Balanced dataset (1:1)

***

## Step 3 — Model Training

LS-GKM parameters:

```
l = 10
k = 6
d = 3
T = 6
```

Model file:

```
models/H3K27me3_VSMC.model.txt
```

***

## Step 4 — Model Performance

| Metric | Value      |
| ------ | ---------- |
| AUROC  | **0.9945** |
| AUPRC  | **0.9946** |

***

## Step 5 — Variant Prioritisation

* Input: cardiovascular GWAS SNPs
* For each variant:
  * Extract **200 bp window**
  * Generate **REF and ALT sequences**
  * Score using LS-GKM

### Delta-score:

```
delta = ALT_score − REF_score
```

Interpretation:

| Delta    | Meaning                            |
| -------- | ---------------------------------- |
| Positive | Gain of H3K27me3 (repression ↑)    |
| Negative | Loss of H3K27me3 (de-repression ↑) |

***

## Step 6 — TF Motif Analysis

* Tool: **FIMO (MEME Suite)**
* Database: **JASPAR CORE**

For each variant:

* Scan REF and ALT sequences
* Identify motif hits near SNP
* Compute:

```
Δ motif score = ALT − REF
```

***

## Step 7 — Figures Generated

* Figure 1 — Peak reproducibility
* Figure 2 — ROC curve
* Figure 3 — Precision-Recall curve
* Figure 4 — Score distribution
* Figure 5 — Variant prioritisation
* Figure 6 — TF motif disruption

***

# Repository Structure

```
├── run_lsgkm_pipeline.sh        # Main training pipeline
├── variant_prioritisation.py    # SNP scoring
├── run_tf_fimo.py              # TF motif analysis
├── generate_figures_and_tf.py  # Figure generation
├── fix_fimo.py                 # Script fix for FIMO
│
├── models/
│   └── H3K27me3_VSMC.model.txt
│
├── variants/
│   ├── snp_ref_windows.fa
│   ├── snp_alt_windows.fa
│   └── variant_prioritisation_results.tsv
│
├── tf_analysis/
│   ├── tf_motif_disruption.tsv
│   └── figures/
│
├── scores/
│   ├── pos_test_scores.txt
│   └── neg_test_scores.txt
│
└── figures/
```

***

# Requirements

## Software

* Python 3.10+
* Bedtools
* MEME Suite (FIMO)
* LS-GKM (gkmtrain, gkmpredict)

## Python packages

```bash
pip install numpy pandas matplotlib seaborn scikit-learn biopython
```

***

# Usage

## 1. Train model

```bash
chmod +x run_lsgkm_pipeline.sh
./run_lsgkm_pipeline.sh
```

***

## 2. Variant prioritisation

```bash
python3 variant_prioritisation.py
```

***

## 3. TF motif analysis

```bash
conda activate meme
python3 run_tf_fimo.py
```

***

## 4. Generate figures

```bash
python3 generate_figures_and_tf.py
```

***

# Example Results

Top variant effects:

| Variant    | Gene       | Delta  | Effect             |
| ---------- | ---------- | ------ | ------------------ |
| rs1260326  | GCKR       | +0.070 | Gain of repression |
| rs1746048  | CXCL12     | -0.066 | Loss of repression |
| rs4977574  | CDKN2B-AS1 | +0.059 | Gain               |
| rs58542926 | TM6SF2     | -0.059 | Loss               |

***

# Biological Interpretation

* Variants can **modulate chromatin repression (H3K27me3)**
* Effects are mediated through:
  * altering TF binding
  * changing motif affinity
* Links GWAS SNPs → **functional regulatory mechanisms in VSMCs**

***

# Methods (for manuscript)

> H3K27me3 ChIP-seq data from VSMCs were processed using nf-core/chipseq aligned to GRCh38. Peaks reproducible in at least 2 of 3 replicates were used as positives. Sequences were centred and resized to 200 bp. An equal number of background sequences were sampled. An LS-GKM model (l=10, k=6, d=3) was trained and evaluated (AUROC=0.995). Variants were prioritised based on delta-scores between reference and alternate alleles.

***

# Author

**Maryamu Usman**  
Demonstrator, University of Edinburgh

***

# Citation

If you use this pipeline:

```
Usman M. (2026)
H3K27me3 LS‑GKM Variant Prioritisation Pipeline
```

***

# Notes

* Works with **Ensembl GRCh38 (no chr prefix)**
* Requires **preprocessed nf-core ChIP-seq output**
* Designed for **epigenomic variant prioritisation**

***

