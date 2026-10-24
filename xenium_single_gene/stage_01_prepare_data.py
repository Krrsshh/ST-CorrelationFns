
# ============================================================
# PIPELINE 5 (XENIUM)
# STAGE 1 — Load Data and Prepare Core-Level Hepatocyte Objects
# ============================================================

import scanpy as sc
import numpy as np


# ============================================================
# Configuration
# ============================================================

INPUT_H5AD = (
    "/Users/krishna/Desktop/LiverFibrosis/Jorge_Mouse_Data/"
    "adata_liver_corrected_LOCALFIX.h5ad"
)

N_TOP_HVGS = 5500

EXPRESSION_LAYER = "lognorm"
COUNTS_LAYER = "counts"

RANDOM_STATE = 42


# ============================================================
# Core-to-sample mapping
# ============================================================
# Each core is assigned to its experimental condition and
# corresponding sample identifier.
#
# Multiple cores may belong to the same sample.
# These mappings are used to organize core-level analyses.

CONTROL_CORES = {
    "core_0":  "Liver_Ctrl_R1.1",
    "core_1":  "Liver_Ctrl_R1.1",
    "core_17": "Liver_Ctrl_R1.2",
    "core_4":  "Liver_Ctrl_R3.1",
    "core_3":  "Liver_Ctrl_R3.1",
    "core_19": "Liver_Ctrl_R3.2",
    "core_2":  "Liver_Ctrl_R3.2",
}

CCL4_CORES = {
    "core_18": "Liver_Ccl4_R2.2",
    "core_6":  "Liver_Ccl4_R2.2",
    "core_8":  "Liver_Ccl4_R4",
}


# ============================================================
# Stage 1 — Prepare analysis inputs
# ============================================================

def prepare_data(
    input_h5ad=INPUT_H5AD,
    n_top_hvgs=N_TOP_HVGS,
    expression_layer=EXPRESSION_LAYER,
    counts_layer=COUNTS_LAYER,
    random_state=RANDOM_STATE,
):
    """
    Load the Xenium AnnData object and prepare hepatocyte-only,
    core-level inputs for the single-gene spatial correlation
    analysis.

    Hepatocyte-specific highly variable genes (HVGs) are selected
    across all retained hepatocytes using raw counts.

    Returns
    -------
    core_data : dict
        Core-level AnnData objects and their condition/sample labels.
    hep_hvg_genes : list
        Gene names selected for downstream analysis.
    """

    # --------------------------------------------------------
    # Reproducibility
    # --------------------------------------------------------
    np.random.seed(random_state)

    # --------------------------------------------------------
    # Load AnnData
    # --------------------------------------------------------
    print("\n======================================")
    print("Loading Xenium AnnData")
    print("======================================")

    adata = sc.read_h5ad(input_h5ad)

    print(adata)

    # --------------------------------------------------------
    # Keep hepatocytes only
    # --------------------------------------------------------
    # All downstream gene selection and spatial correlation
    # calculations are restricted to annotated hepatocytes.

    required_obs_columns = ["cell_type", "core_id_raw"]

    missing_columns = [
        col for col in required_obs_columns
        if col not in adata.obs.columns
    ]

    if missing_columns:
        raise KeyError(
            f"Missing required AnnData.obs columns: {missing_columns}"
        )

    required_layers = [expression_layer, counts_layer]

    missing_layers = [
        layer for layer in required_layers
        if layer not in adata.layers
    ]

    if missing_layers:
        raise KeyError(
            f"Missing required AnnData layers: {missing_layers}"
        )

    adata_hep = adata[
        adata.obs["cell_type"] == "Hepatocyte"
    ].copy()

    print("\n======================================")
    print("HEP-only object")
    print("======================================")
    print(adata_hep)

    print("\nHEP cells per core:")
    print(adata_hep.obs["core_id_raw"].value_counts())

    # --------------------------------------------------------
    # Compute hepatocyte-specific highly variable genes
    # --------------------------------------------------------
    # HVGs are selected jointly across all retained hepatocytes,
    # rather than independently within each core.
    #
    # Raw counts are used for feature selection. The resulting
    # gene list is shared across all cores and both conditions.

    print("\n======================================")
    print("Computing HEP-specific HVGs")
    print("======================================")

    adata_hep.X = adata_hep.layers[counts_layer].copy()

    try:
        sc.pp.highly_variable_genes(
            adata_hep,
            n_top_genes=n_top_hvgs,
            flavor="seurat_v3",
            layer=counts_layer,
        )

    except Exception as e:
        # Preserve the original fallback behavior if Seurat-v3
        # HVG selection fails in the current environment.
        print("\nseurat_v3 failed")
        print("Using cell_ranger instead")
        print(e)

        sc.pp.highly_variable_genes(
            adata_hep,
            n_top_genes=n_top_hvgs,
            flavor="cell_ranger",
        )

    hep_hvg_genes = (
        adata_hep.var_names[
            adata_hep.var["highly_variable"]
        ].tolist()
    )

    print(f"\nSelected HEP HVGs: {len(hep_hvg_genes)}")

    # Restore log-normalized expression for downstream
    # spatial correlation calculations.
    adata_hep.X = adata_hep.layers[expression_layer].copy()

    # --------------------------------------------------------
    # Build core-level data dictionary
    # --------------------------------------------------------
    # Each entry contains:
    #   condition : "Control" or "CCl4"
    #   sample_id : corresponding sample identifier
    #   adata     : hepatocyte-only AnnData object for that core
    #
    # Core eligibility based on cell count and spatial-pair
    # availability is evaluated later in Stage 2.

    core_data = {}

    for condition, core_dict in [
        ("Control", CONTROL_CORES),
        ("CCl4", CCL4_CORES),
    ]:
        for core_id, sample_id in core_dict.items():

            sub = adata_hep[
                adata_hep.obs["core_id_raw"] == core_id
            ].copy()

            print(f"\nLoading {condition}: {core_id}")
            print(f"Sample     : {sample_id}")
            print(f"HEP cells  : {sub.n_obs}")

            core_data[core_id] = {
                "condition": condition,
                "sample_id": sample_id,
                "adata": sub,
            }

    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------
    print("\n======================================")
    print("Stage 1 Complete")
    print("======================================")

    print(f"Control cores  : {len(CONTROL_CORES)}")
    print(f"CCl4 cores     : {len(CCL4_CORES)}")
    print(f"Total cores    : {len(core_data)}")
    print(f"Analysis genes : {len(hep_hvg_genes)}")

    return core_data, hep_hvg_genes


# ============================================================
# Allow Stage 1 to run independently for testing
# ============================================================

if __name__ == "__main__":
    core_data, hep_hvg_genes = prepare_data()
