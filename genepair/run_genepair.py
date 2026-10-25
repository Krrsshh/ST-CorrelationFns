
"""
Driver script for the Xenium gene-pair spatial correlation pipeline.

Runs:
    Stage 1: Prepare data
    Stage 2: Compute spatial correlations and fit correlation lengths
    Stage 3: Summarize condition-level effects
"""

from genepair.stage_01_prepare_data import prepare_data
from genepair.stage_02_compute_spatial_correlations import (
    compute_spatial_correlations,
)
from genepair.stage_03_summarize_effects import summarize_effects


def main():
    # Stage 1: Prepare data
    core_data, genes_use, core_summary_df = prepare_data()

    # Stage 2: Compute gene-pair spatial correlations
    xi_batch_pair_df, pair_qc_df, core_paircount_results = (
        compute_spatial_correlations(core_data, genes_use)
    )

    # Stage 3: Summarize condition-level effects
    effect_pair_df = summarize_effects(xi_batch_pair_df)

    print("\nGene-pair pipeline completed successfully.")
    print(f"Genes analyzed: {len(genes_use)}")
    print(f"Gene-pair correlation records: {len(xi_batch_pair_df)}")
    print(f"Final gene-pair effect records: {len(effect_pair_df)}")

    return (
        core_summary_df,
        xi_batch_pair_df,
        pair_qc_df,
        core_paircount_results,
        effect_pair_df,
    )


if __name__ == "__main__":
    main()
