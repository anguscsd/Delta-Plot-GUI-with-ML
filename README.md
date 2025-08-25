# Repeat Expansion Analysis GUI

## Overview

This repository contains a Python-based graphical user interface (GUI) designed to streamline the analysis of CAG repeat expansions, with a focus on Huntington’s disease datasets. The application integrates the Repeat Detector (RD) workflow within a user-friendly environment, enabling users to generate histograms, calculate instability metrics, and produce publication-ready delta plots without requiring command-line interaction or coding experience.

The tool is implemented in **Python 3.13** using **wxPython** for the interface, **Docker** for executing the Repeat Detector workflow, and **pandas**, **numpy**, and **matplotlib** for data handling, statistical calculations, and visualisation.

The GUI provides a fully integrated pipeline that takes raw FASTA files as input and produces interpretable outputs, including histogram files, instability indices, and delta plots, with options for dataset comparison and error bar visualisation.

---

## Key Features

### Integrated Pipeline
The GUI integrates the Repeat Detector workflow and automates the generation of histogram files, instability metrics, and visualisations, providing a seamless process from raw data to results.

### Two Main Analysis Modules

#### 1. Generate Histograms
This tab provides controls for running the Repeat Detector workflow:
* Select the Docker `.tar` file containing the RD environment.
* Specify the input directory containing FASTA files.
* Choose analysis profiles: **restrictive** or **permissive**.
* Optionally enable calculation of instability metrics, including Instability Index (II), Expansion Index (EI), and Contraction Index (CI).
* Set an output directory for histogram files.
* Automatically sorts results into `day0` and `day42` subdirectories based on filename conventions.
* Displays real-time status updates while RD is running.

#### 2. Analyse Histograms
This tab enables statistical analysis and visualisation of histogram outputs:
* Load histogram files for **day 0** and **day 42** for one or two datasets.
* Supports both **tab-separated** and **comma-separated** file formats.
* Normalises read counts using the equation:  
  **normalised_reads = (reads / Σ reads) × 100**
* Calculates delta values:  
  **Δ = Day42_mean − Day0_mean**
* Optionally displays **error bars** between replicates or datasets:
  * Standard Deviation (SD)
  * Standard Error of the Mean (SEM)
* Generates delta plots with customisable bin size, plot title, and colour selection.
* Exports processed datasets and visualisations as CSV and PDF files.
* Embeds plots directly within the GUI and supports external PDF viewing.

---

## Installation

### Prerequisites
Before running the application, ensure you have the following installed:
* **Python 3.13 or later**
* **Docker** (for Repeat Detector execution)
* Required Python libraries:
  ```bash
  pip install wxPython pandas numpy matplotlib

# R Script fot ML
**This current version of the Machine Learning is written and run in R**

The pipeline supports two clustering approaches:
1. **Gaussian Mixture Models (GMM)** *(recommended)*  
2. **K-means Clustering** *(alternative method)*

It generates a comprehensive PDF report containing dataset-specific results, bin-level predictions, and **Receiver Operating Characteristic (ROC)** curves for model evaluation, alongside structured CSV outputs for further analysis.

---

## Input Data Format
The input CSV files should follow this structure:
- **First column**: Repeat length (numeric, ascending or mixed order).
- **Second column**: Δ values (Day42 − Day0).  
- File names must be unique and meaningful, as they are used in result summaries.

---

## Workflow Summary

### 1. Data Preprocessing
- The pipeline reads raw CSV files and bins repeat lengths into intervals defined by `bin_width` (default = 1).
- Within each bin, the mean Δ value is calculated.
- Non-numeric and missing values are excluded.
- Feature scaling is applied via **z-score normalization** using training dataset statistics only.

### 2. Model Training
#### Gaussian Mixture Model (GMM)
- Implemented via the **mclust** package.
- Models the distribution of Δ values as a **mixture of Gaussian components**.
- Supports **soft clustering**, assigning each bin a **probability** of belonging to each cluster.
- Automatically determines the optimal number of clusters (K) when `auto_select_K = TRUE` using the **Bayesian Information Criterion (BIC)**.
- Uses the `tau` parameter to classify bins:
  - Δ > +τ → **Expansion**
  - Δ < −τ → **Contraction**
  - |Δ| ≤ τ → **Neutral**

#### K-means Clustering
- Uses the **stats::kmeans** implementation.
- Clusters are assigned based on minimizing **Euclidean distance** between bins and their nearest cluster centroid.
- Produces **hard clustering** (each bin is placed into a single cluster, no probability weights).
- Less flexible than GMM when cluster shapes are overlapping or uneven but included for comparison.

### 3. Model Evaluation
- Predictions are generated for each bin in the **test datasets** using the trained model.
- **ROC curves** are calculated for both **expansion** and **contraction** detection using the **pROC** package.
- For each test dataset:
  - Two ROC curves are generated (expansion and contraction).
  - **AUC (Area Under Curve)** values are computed with 95% confidence intervals.
- An **overall ROC curve** summarizing performance across all test datasets is also produced.

### 4. Outputs
The pipeline produces the following outputs in the `unsup_model_out/` directory:
- **PDF Report**  
  - Multi-panel plots per dataset (observed Δ, predicted Δ, ROC curves).
  - Overall ROC summary figure.
- **CSV Files**  
  - `train_assignments.csv` → Cluster assignments for training data.
  - `cluster_centers.csv` → Cluster centers in Δ space.
  - `test_metrics.csv` → AUC values per dataset.
  - `test_bin_predictions.csv` → Predicted Δ values per bin.

---

## Key Parameters
| Parameter        | Description                                    | Default |
|------------------|-----------------------------------------------|---------|
| `method`        | Clustering algorithm: `"gmm"` or `"kmeans"`    | `"gmm"` |
| `bin_width`     | Repeat-length bin size                        | `1`     |
| `tau`           | Threshold for expansion/contraction labeling  | `0.0075` |
| `auto_select_K` | Auto-select number of clusters (GMM only)     | `TRUE`  |
| `smooth_k`      | Rolling mean smoothing window for plots       | `5`     |

---

## Dependencies
The following R packages are required:
- **tidyverse** → Data handling
- **ggplot2** → Visualization
- **zoo** → Smoothing observed Δ signals
- **pROC** → ROC curve computation and AUC metrics
- **patchwork** → Combining multi-panel figures
- **scales** → Formatting axes and labels
- **readr** → Fast CSV reading/writing
- **mclust** → Gaussian Mixture Models

Install automatically on first run if not present.

---

## Recommended Use
- **Default**: Use `method = "gmm"` for more robust clustering of overlapping Δ distributions.
- **When to use K-means**: If clusters are well-separated, spherical, and non-overlapping.
- For small datasets (<10 files), results may be less stable; performance improves as more datasets are added.

---

## Summary
This pipeline provides a reproducible framework for detecting **repeat expansion and contraction patterns** using **unsupervised machine learning**. By combining GMM-based probabilistic modeling, bin-level Δ analysis, and ROC-based evaluation, it enables sensitive detection of subtle repeat instability trends across datasets. The integration of both GMM and K-means offers flexibility, though GMM is recommended for datasets where biological variation produces overlapping distributions.

---

**Citation**
If using this pipeline in publications, cite:
> "Gaussian Mixture Models and unsupervised learning were implemented using the `mclust` R package (Scrucca et al., 2016), and ROC analysis performed with the `pROC` package (Robin et al., 2011)."
