
# ============================================================
# XENIUM GENE–GENE CROSS-CORRELATION PIPELINE
# STAGE 1 — Load Core Objects and Select Analysis Genes
#
# Healthy HEP vs Diseased HEP
# ============================================================

import scanpy as sc
import pandas as pd
import numpy as np


# ============================================================
# Paths
# ============================================================

INPUT_H5AD = (
    "./adata_liver_corrected_LOCALFIX.h5ad"
)

ORTHOLOGUE_FILE = (
    "./MERFISH_human_mouse_orthologues.xlsx"
)


# ============================================================
# Core configuration
# ============================================================

healthy_cores = {

    "core_0":  "Liver_Ctrl_R1.1",
    "core_1":  "Liver_Ctrl_R1.1",
    "core_17": "Liver_Ctrl_R1.2",
    "core_4":  "Liver_Ctrl_R3.1",
    "core_3":  "Liver_Ctrl_R3.1",
    "core_19": "Liver_Ctrl_R3.2",
    "core_2":  "Liver_Ctrl_R3.2"
}


diseased_cores = {

    "core_18": "Liver_Ccl4_R2.2",
    "core_6":  "Liver_Ccl4_R2.2",
    "core_8":  "Liver_Ccl4_R4"
}


# ============================================================
# Parameters
# ============================================================

expression_layer = "lognorm"

USE_MERFISH_ORTHOLOGUES = False

random_state = 42
np.random.seed(random_state)


# ============================================================
# Main Stage 1 function
# ============================================================

def prepare_data():

    # ========================================================
    # STEP 1 — Load Xenium object
    # ========================================================

    print("\n==============================================")
    print("STAGE 1 — Loading Xenium HEP object")
    print("==============================================")

    adata = sc.read_h5ad(INPUT_H5AD)

    print(adata)


    # ========================================================
    # Keep hepatocytes only
    # ========================================================

    adata_hep = adata[
        adata.obs["cell_type"] == "Hepatocyte"
    ].copy()

    adata_hep.X = adata_hep.layers[
        expression_layer
    ].copy()

    print("\n==============================================")
    print("HEP object")
    print("==============================================")

    print(adata_hep)

    print("\nHEP cells per core:")

    print(
        adata_hep.obs["core_id_raw"]
        .value_counts()
    )


    # ========================================================
    # STEP 2 — Select analysis genes
    # ========================================================

    print("\n==============================================")
    print("Selecting analysis genes")
    print("==============================================")

    if USE_MERFISH_ORTHOLOGUES:

        # ----------------------------------------------------
        # Load MERFISH mouse orthologue table
        # ----------------------------------------------------

        orth_df = pd.read_excel(ORTHOLOGUE_FILE)

        # ----------------------------------------------------
        # Extract unique mouse orthologues
        # ----------------------------------------------------

        mouse_genes = []

        for entry in orth_df["Mouse orthologue"].dropna():

            genes = [
                g.strip()
                for g in str(entry).split(",")
            ]

            mouse_genes.extend(genes)

        mouse_genes = sorted(set(mouse_genes))

        # ----------------------------------------------------
        # Keep only genes present in Xenium
        # ----------------------------------------------------

        genes_use = sorted(
            list(
                set(mouse_genes)
                &
                set(adata_hep.var_names)
            )
        )

        print(f"Unique mouse orthologues: {len(mouse_genes)}")
        print(f"Mouse orthologues present in Xenium: {len(genes_use)}")

    else:

        genes_use = (
            adata_hep.var_names
            .tolist()
        )

        print(f"Using all Xenium genes: {len(genes_use)}")


    print("\n==============================================")
    print("Selected analysis genes")
    print("==============================================")

    print(f"Genes selected: {len(genes_use)}")

    print(
        f"Possible unordered gene pairs: "
        f"{len(genes_use) * (len(genes_use)-1)//2:,}"
    )


    # --------------------------------------------------------
    # Keep only selected genes
    # --------------------------------------------------------

    adata_hep = adata_hep[
        :,
        genes_use
    ].copy()

    adata_hep.X = adata_hep.layers[
        expression_layer
    ].copy()


    # ========================================================
    # STEP 3 — Load all cores
    # ========================================================

    core_data = {}

    core_summary_rows = []

    print("\n==============================================")
    print("STAGE 1 — Loading HEP cores")
    print("==============================================")


    for condition, core_dict in [

        ("Healthy", healthy_cores),

        ("Diseased", diseased_cores)

    ]:

        for core_id, sample_id in core_dict.items():

            print("\n----------------------------------------------")
            print(
                f"Loading {condition}: {core_id}"
            )
            print("----------------------------------------------")

            adata_core = adata_hep[
                adata_hep.obs["core_id_raw"] == core_id
            ].copy()

            core_data[core_id] = {

                "condition": condition,

                "sample_id": sample_id,

                "adata": adata_core,

                "n_hepatocytes": adata_core.n_obs

            }

            core_summary_rows.append({

                "Core": core_id,

                "Sample": sample_id,

                "Condition": condition,

                "n_hepatocytes": adata_core.n_obs,

                "n_expression_genes": adata_core.n_vars

            })

            print(f"HEP cells: {adata_core.n_obs:,}")
            print(f"Expression genes: {adata_core.n_vars:,}")


    # ========================================================
    # Core summary
    # ========================================================

    core_summary_df = pd.DataFrame(
        core_summary_rows
    )

    print("\n==============================================")
    print("Core loading summary")
    print("==============================================")

    print(core_summary_df.to_string(index=False))


    # ========================================================
    # STEP 4 — Confirm all selected genes exist in every core
    # ========================================================

    for core_id, info in core_data.items():

        missing_genes = [
            gene
            for gene in genes_use
            if gene not in info["adata"].var_names
        ]

        if missing_genes:

            raise ValueError(
                f"{core_id} is missing selected genes: "
                f"{missing_genes}"
            )


    # ========================================================
    # Final Stage 1 summary
    # ========================================================

    n_healthy_cells = sum(

        info["n_hepatocytes"]

        for info in core_data.values()

        if info["condition"] == "Healthy"
    )


    n_diseased_cells = sum(

        info["n_hepatocytes"]

        for info in core_data.values()

        if info["condition"] == "Diseased"
    )


    print("\n==============================================")
    print("STAGE 1 COMPLETE")
    print("==============================================")

    print(f"Healthy cores: {len(healthy_cores)}")
    print(f"Diseased cores: {len(diseased_cores)}")
    print(f"Total cores: {len(core_data)}")

    print(f"Healthy HEP cells: {n_healthy_cells:,}")
    print(f"Diseased HEP cells: {n_diseased_cells:,}")

    print(f"Analysis genes: {len(genes_use):,}")

    print(
        f"Gene pairs per core: "
        f"{len(genes_use) * (len(genes_use)-1)//2:,}"
    )

    print(
        f"Total core–gene-pair combinations: "
        f"{len(core_data) * len(genes_use) * (len(genes_use)-1)//2:,}"
    )

    # ========================================================
    # Return objects required by subsequent stages
    # ========================================================

    return core_data, genes_use, core_summary_df
