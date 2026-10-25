# Spatial Correlation Functions in Transcriptomics

## 1. Project Overview

This repository contains computational pipelines for implementing spatial gene-expression correlation functions in spatial transcriptomics data. The approach characterizes how gene-expression correlations vary with the distance between cells, providing measures of spatial correlation length $(\xi)$ and susceptibility $(\chi)$. The repository currently includes a single-gene analysis pipeline and a gene-pair analysis pipeline for hepatocytes in mouse liver Xenium data, comparing Control and CCl4-treated samples.

## 2. Repository Contents

The repository contains two analysis pipelines: a single-gene pipeline implemented in `xenium_single_gene/` and a gene-pair pipeline implemented in `genepair/`.

### Single-Gene Pipeline

```text
xenium_single_gene/
├── __init__.py
├── run_xenium_single_gene.py
├── stage_01_prepare_data.py
├── stage_02_compute_spatial_correlations.py
└── stage_03_summarize_effects.py
```

- **`run_xenium_single_gene.py`** — Runs the three analysis stages and exports the final effect table.
- **`stage_01_prepare_data.py`** — Loads the AnnData dataset, selects hepatocytes, identifies highly variable genes (HVGs), and prepares the configured cores.
- **`stage_02_compute_spatial_correlations.py`** — Computes core-level gene-expression correlation functions, estimates $(\xi)$, and calculates $(\chi)$.
- **`stage_03_summarize_effects.py`** — Summarizes core-level results by condition and calculates Control-versus-CCl4 effect metrics.
- **`__init__.py`** — Initializes the Python package.

### Gene-Pair Pipeline

```text
genepair/
├── __init__.py
├── run_genepair.py
├── stage_01_prepare_data.py
├── stage_02_compute_spatial_correlations.py
└── stage_03_summarize_effects.py
```

- **`run_genepair.py`** — Runs the three analysis stages and returns the final gene-pair effect table and intermediate results.
- **`stage_01_prepare_data.py`** — Loads the AnnData dataset, selects hepatocytes, and prepares the configured cores and genes for gene-pair analysis. By default, all genes are used.
- **`stage_02_compute_spatial_correlations.py`** — Computes core-level gene-pair expression correlation functions, estimates $(\xi)$, and calculates $(\chi)$.
- **`stage_03_summarize_effects.py`** — Summarizes gene-pair results by condition and calculates Control-versus-CCl4 effect metrics.
- **`__init__.py`** — Initializes the Python package.

The `outputs/` directory is used to store generated result files.

## 3. Method Overview

The analysis is performed separately for each gene in each core.

1. **Spatial pairs:** A KD-tree is used to identify unique, unordered pairs of hepatocytes separated by no more than 100 µm.
2. **Correlation function:**
   - **Single-gene pipeline:** Computes spatial autocorrelation for each gene using mean-centered expression products normalized by gene-expression variance.
   - **Gene-pair pipeline:** Computes spatial cross-correlation between pairs of genes using mean-centered expression products normalized by the product of their gene-expression standard deviations.
3. **Correlation length $(\xi)$:** An exponential decay model, $(C(r)=A\exp(-r/\xi))$, is fitted to the correlation function using distance-bin centers at or beyond 30 µm. The fitted $(\xi)$ characterizes the spatial decay length of the correlation.
4. **Susceptibility $(\chi)$:** The correlation function is integrated over the same distance tail using radial weighting, $(2\pi r\Delta r\)$.
5. **Condition-level effects:** Core-level results are summarized separately for Control and CCl4, and effect metrics are calculated to compare the conditions.

The distance-bin width is 5 µm. The analysis uses the `lognorm` expression layer by default.

## 4. Filtering and Quality Control

Filtering occurs at three levels: core selection during data preparation, core-gene or core-gene-pair evaluation during spatial correlation analysis, and gene or gene-pair selection during condition-level summarization.

### Stage 1 — Gene and Core Selection

- Only hepatocytes (`cell_type == "Hepatocyte"`) are retained.
- Single-gene pipeline: Analysis is restricted to 5,500 highly variable genes selected from hepatocyte raw counts using Scanpy's seurat_v3 method, with cell_ranger as a fallback.
- Gene-pair pipeline: All genes are used by default. HVG selection could be implemented similarly if required.
- Ten configured cores are considered: seven Control cores and three CCl4 cores.
- Cores containing fewer than 100 hepatocytes are skipped during Stage 2.
- Cores without valid spatial pairs within 100 µm are skipped.

### Stage 2 — Core-Level Gene-Wise $(\xi)$ and $(\chi)$

The following QC checks are applied to individual genes or gene pairs, as appropriate:

- **Gene variance:** Expression variance must be finite and greater than zero.
- **Finite correlation curve:** At least one finite correlation value must be available.
- **Minimum correlation magnitude:** `min_abs_corr = 0.0`; this setting does not impose an effective additional filter. Only retained for the user if user wants a more stricter result filtering.
- **Minimum tail points:** At least five finite correlation values must be available at distance-bin centers of 30 µm or greater.
- **Curve variability:** The standard deviation of the full correlation curve must be at least $(10^{-4})$.

The exponential fit is accepted only if numerical fitting succeeds, $(\xi>0)$, $(\xi<500)$ µm, and $(R^2>0.1)$. The fit is unconstrained, the positive-tail requirement is disabled, and the estimated uncertainty in $(\xi)$ is not used as an acceptance criterion.

**Retention of $(\chi)$:** The integrated correlation is calculated over the same distance tail. If fitting fails or the fitted $(\xi)$ is rejected, the result can still be stored with $(\xi)$ recorded as missing (`NaN`). However, combinations excluded by the earlier variance, correlation-curve, tail-point, or curve-variability checks do not contribute to either $(\xi)$ or $(\chi)$ summaries.

### Stage 3 — Condition-Level Effect Table

Condition-level summaries use the following criteria:

- **Valid $(\xi)$ estimates:** Only rows with accepted, non-missing $(\xi)$ contribute to condition-level $(\xi)$ means, standard deviations, counts, and mean $\(R^2)$.
- **$(\chi)$ summaries:** $(\chi)$ is summarized independently of $(\xi)$-fit acceptance, using stored rows with finite $(\chi)$.
- **Condition coverage:** Both Control and CCl4 must have usable summaries for the required metrics.
- **Positive mean correlation length:** Mean $(\xi)$ must be positive in both conditions.
- **Minimum replication:** At least three valid core-level $(\xi)$ estimates are required in each condition.
- **Finite effect metrics:** Only rows with finite values for all required effect metrics are retained.

The effect metrics differ between the pipelines:

- Single-gene pipeline: Calculates $(Z_\xi)$, $(Z_\chi)$, $(\log(\xi_{\mathrm{Control}}/\xi_{\mathrm{CCl4}}))$, and $(\Delta\chi)$.
- Gene-pair pipeline: Calculates $(Z_\xi)$, $(Z_\chi)$, $(Z_{|\chi|})$, $(\log(\xi_{\mathrm{Control}}/\xi_{\mathrm{CCl4}}))$, $(\Delta\chi)$, and $(\Delta|\chi|)$.

No additional $(R^2)$ threshold is applied in Stage 3. Mean $(R^2)$ values are reported as descriptive quality-control summaries of accepted fits.

**Implementation note:** Although $(\chi)$ does not depend on whether the $(\xi)$ fit is accepted, its inclusion still depends on the gene-core combination passing the preceding Stage 2 checks and reaching the result-storage step.

## 5. Input Data and Prerequisites

### Input data

The pipeline expects an AnnData (`.h5ad`) object containing:

- **`adata.obs["cell_type"]`** — Cell-type annotations used to select hepatocytes.
- **`adata.obs["core_id_raw"]`** — Core identifiers used to assign cells to the configured Control and CCl4 cores.
- **`adata.layers["counts"]`** — Raw gene counts used for highly variable gene selection.
- **`adata.layers["lognorm"]`** — Log-normalized expression values used for spatial correlation analysis.
- **`adata.obsm["spatial"]`** — Two-dimensional spatial coordinates for the cells, expressed in micrometers.

The gene identifiers in the expression matrix must be compatible with the gene list produced during data preparation. The configured core identifiers and their condition assignments are defined in `stage_01_prepare_data.py`.

The input dataset is not included in this repository. Users must provide their own appropriately processed dataset and configure the input path in `stage_01_prepare_data.py`.

### Software requirements

The pipeline was developed and tested with the following Python package versions:

- Python 3
- Scanpy 1.11.5
- AnnData 0.11.4
- NumPy 1.26.3
- pandas 2.2.0
- SciPy 1.12.0
- scikit-learn 1.4.1.post1

These dependencies are listed in `requirements.txt`.

## 6. How to Run

### Step 1 — Set up the environment

Clone or download the repository, then navigate to its root directory.

Create and activate a Python environment if needed, and install the dependencies:

```bash
python -m pip install -r requirements.txt
```

### Step 2 — Configure the input data

Configure the input .h5ad path in the relevant stage_01_prepare_data.py file:
- Single-gene pipeline: xenium_single_gene/stage_01_prepare_data.py
- Gene-pair pipeline: genepair/stage_01_prepare_data.py
Ensure that the required metadata columns, expression layers, and spatial coordinates described in Section 5 are available.

### Step 3 — Run the pipeline

From the repository root, execute the command for the desired analysis.

Single-gene pipeline:
```bash
python -m xenium_single_gene.run_xenium_single_gene
```
Gene-pair pipeline:
```bash
python -m genepair.run_genepair
```

Each driver runs the three stages in order:
1. Data preparation and gene selection.
2. Core-level spatial correlation analysis.
3. Condition-level effect summarization.
The single-gene pipeline selects highly variable genes (HVGs), whereas the gene-pair pipeline uses all genes by default.

### Step 4 — Inspect the output

The final effect tables are written to:

```text
outputs/Corrected_Xenium_SingleGene_EffectTable_y0relaxed.csv
outputs/Corrected_Xenium_GenePair_EffectTable.csv
```

The tables contains condition-level summaries of $(\xi)$ and $(\chi)$, sample counts, mean fit-quality values, and the following effect metrics:

- `Z_xi` — Standardized difference in mean correlation length.
- `Z_chi` — Standardized difference in mean integrated correlation.
- `log_xi_ratio` — Logarithm of the ratio of mean correlation lengths, Control relative to CCl4.
- `delta_chi` — Difference in mean integrated correlation, Control minus CCl4.

Intermediate core-level results and quality-control information are returned by the analysis stages. The driver currently exports the final condition-level effect table.

## 7. Reproducibility and Citation
The single-gene pipeline uses Scanpy for highly variable gene (HVG) selection and sets the HVG selection random state to 42. Relevant software dependencies are recorded in `requirements.txt`. Spatial analysis settings, core assignments, and filtering thresholds are defined in the corresponding stage scripts.

For reproducibility, use a compatible input AnnData object and retain the documented metadata, expression layers, spatial coordinates, and core assignments. Results may depend on the input dataset and preprocessing.

**Manuscript:** *A correlation-function approach to spatial transcriptomics reveals zonation changes in liver disease.* Full citation details will be added when available.

**License:** A license has not yet been specified. Terms for reuse and redistribution should be confirmed with the research group before public release.

