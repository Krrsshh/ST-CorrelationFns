
# ============================================================
# PIPELINE 5 (XENIUM)
# Driver — Single-Gene Spatial Correlation Analysis
# ============================================================

from pathlib import Path

from xenium_single_gene.stage_01_prepare_data import prepare_data
from xenium_single_gene.stage_02_compute_spatial_correlations import (
    compute_spatial_correlations,
)
from xenium_single_gene.stage_03_summarize_effects import summarize_effects


# ============================================================
# Configuration
# ============================================================

EXPRESSION_LAYER = "lognorm"

OUTPUT_DIR = Path("outputs")

OUTPUT_EFFECT_TABLE = (
    OUTPUT_DIR
    / "Corrected_Xenium_SingleGene_EffectTable_y0relaxed.csv"
)


# ============================================================
# Stage 1 — Prepare data
# ============================================================

print("\n" + "=" * 60)
print("STAGE 1 — DATA PREPARATION")
print("=" * 60)

core_data, hep_hvg_genes = prepare_data()

print(f"\nPrepared cores: {len(core_data)}")
print(f"Hepatocyte HVGs: {len(hep_hvg_genes)}")


# ============================================================
# Stage 2 — Compute core-level spatial correlations
# ============================================================

print("\n" + "=" * 60)
print("STAGE 2 — SPATIAL CORRELATION ANALYSIS")
print("=" * 60)

xi_core_gene_df, qc = compute_spatial_correlations(
    core_data=core_data,
    hep_hvg_genes=hep_hvg_genes,
    expression_layer=EXPRESSION_LAYER,
)

print(f"\nCore-level result rows: {len(xi_core_gene_df)}")


# ============================================================
# Stage 3 — Summarize condition-level effects
# ============================================================

print("\n" + "=" * 60)
print("STAGE 3 — CONDITION-LEVEL EFFECT TABLE")
print("=" * 60)

effect_df = summarize_effects(xi_core_gene_df)


# ============================================================
# Export final effect table
# ============================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

effect_df.to_csv(
    OUTPUT_EFFECT_TABLE,
    index=False,
)

print("\n" + "=" * 60)
print("PIPELINE COMPLETE")
print("=" * 60)

print(f"Prepared cores       : {len(core_data)}")
print(f"Hepatocyte HVGs      : {len(hep_hvg_genes)}")
print(f"Core-level result rows: {len(xi_core_gene_df)}")
print(f"Final genes exported : {len(effect_df)}")
print(f"Effect table saved to: {OUTPUT_EFFECT_TABLE.resolve()}")
