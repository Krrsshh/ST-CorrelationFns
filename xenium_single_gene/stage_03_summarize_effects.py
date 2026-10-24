
# ============================================================
# PIPELINE 5 (XENIUM)
# STAGE 3 — Condition-level Effect Table
# Control vs CCl4
# ============================================================

import numpy as np
import pandas as pd


# ============================================================
# Main summarization
# ============================================================

def summarize_effects(xi_core_gene_df):
    """
    Summarize core-level ξ and χ estimates by condition and
    calculate condition-level gene effect metrics.

    Parameters
    ----------
    xi_core_gene_df : pandas.DataFrame
        Output from Stage 2, containing core-level gene-wise
        ξ, χ, and fit statistics.

    Returns
    -------
    effect_df : pandas.DataFrame
        Final condition-level effect table.
    """

    required_input_columns = [
        "Gene",
        "Condition",
        "xi",
        "chi",
        "r2",
    ]

    missing_columns = [
        col for col in required_input_columns
        if col not in xi_core_gene_df.columns
    ]

    if missing_columns:
        raise KeyError(
            f"Stage 2 results are missing required columns: "
            f"{missing_columns}"
        )

    # ========================================================
    # A. ξ summary from valid ξ fits
    # ========================================================

    xi_valid_df = (
        xi_core_gene_df
        .dropna(subset=["xi"])
        .copy()
    )

    xi_summary = (
        xi_valid_df
        .groupby(["Gene", "Condition"])
        .agg(
            xi_mean=("xi", "mean"),
            xi_sd=("xi", "std"),
            xi_n=("xi", "count"),
        )
        .reset_index()
    )

    xi_pivot = xi_summary.pivot(
        index="Gene",
        columns="Condition",
        values=["xi_mean", "xi_sd", "xi_n"],
    )

    xi_pivot.columns = [
        f"{metric}_{condition}"
        for metric, condition in xi_pivot.columns
    ]

    xi_pivot = xi_pivot.reset_index()

    # ========================================================
    # A2. Mean R² summary (QC annotation only)
    # ========================================================

    r2_summary = (
        xi_valid_df
        .groupby(["Gene", "Condition"])
        .agg(r2_mean=("r2", "mean"))
        .reset_index()
    )

    r2_pivot = r2_summary.pivot(
        index="Gene",
        columns="Condition",
        values="r2_mean",
    )

    r2_pivot.columns = [
        f"r2_mean_{condition}"
        for condition in r2_pivot.columns
    ]

    r2_pivot = r2_pivot.reset_index()

    print("\n========================================")
    print("Stage 3 QC")
    print("========================================")
    print(f"Genes with at least one valid ξ fit: {len(xi_pivot)}")

    # ========================================================
    # B. χ summary independent of ξ-fit success
    # ========================================================

    chi_valid_df = (
        xi_core_gene_df
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["chi"])
        .copy()
    )

    chi_summary = (
        chi_valid_df
        .groupby(["Gene", "Condition"])
        .agg(
            chi_mean=("chi", "mean"),
            chi_sd=("chi", "std"),
            chi_n=("chi", "count"),
        )
        .reset_index()
    )

    chi_pivot = chi_summary.pivot(
        index="Gene",
        columns="Condition",
        values=["chi_mean", "chi_sd", "chi_n"],
    )

    chi_pivot.columns = [
        f"{metric}_{condition}"
        for metric, condition in chi_pivot.columns
    ]

    chi_pivot = chi_pivot.reset_index()

    print(f"Genes with valid χ values: {len(chi_pivot)}")

    # ========================================================
    # C. Merge ξ and χ summaries
    # ========================================================

    effect_df = (
        xi_pivot
        .merge(r2_pivot, on="Gene", how="left")
        .merge(chi_pivot, on="Gene", how="inner")
    )

    print(f"After ξ–χ merge: {len(effect_df)} genes")

    required_columns = [
        "xi_mean_Control",
        "xi_mean_CCl4",
        "xi_sd_Control",
        "xi_sd_CCl4",
        "xi_n_Control",
        "xi_n_CCl4",
        "chi_mean_Control",
        "chi_mean_CCl4",
    ]

    effect_df = (
        effect_df
        .dropna(subset=required_columns)
        .copy()
    )

    print(
        f"After requiring both conditions: "
        f"{len(effect_df)} genes"
    )

    # --------------------------------------------------------
    # Require positive condition-level mean ξ values
    # --------------------------------------------------------

    effect_df = effect_df[
        (effect_df["xi_mean_Control"] > 0)
        & (effect_df["xi_mean_CCl4"] > 0)
    ].copy()

    print(
        f"After positive mean ξ filter: "
        f"{len(effect_df)} genes"
    )

    # --------------------------------------------------------
    # Require at least 3 valid cores in each condition
    # --------------------------------------------------------

    effect_df = effect_df[
        (effect_df["xi_n_Control"] >= 3)
        & (effect_df["xi_n_CCl4"] >= 3)
    ].copy()

    print(
        f"After minimum core count filter: "
        f"{len(effect_df)} genes"
    )

    # ========================================================
    # D. Calculate paper-style effect metrics
    # ========================================================

    # Standardized ξ difference:
    # Z_xi = (mean_Control - mean_CCl4)
    #        / sqrt(sd_Control² + sd_CCl4²)

    z_xi_denominator = np.sqrt(
        effect_df["xi_sd_Control"] ** 2
        + effect_df["xi_sd_CCl4"] ** 2
    )

    effect_df["Z_xi"] = (
        effect_df["xi_mean_Control"]
        - effect_df["xi_mean_CCl4"]
    ) / z_xi_denominator

    # Scale-independent ξ effect magnitude.
    effect_df["log_xi_ratio"] = np.log(
        effect_df["xi_mean_Control"]
        / effect_df["xi_mean_CCl4"]
    )

    # Model-free susceptibility difference.
    effect_df["delta_chi"] = (
        effect_df["chi_mean_Control"]
        - effect_df["chi_mean_CCl4"]
    )

    # Standardized susceptibility difference.
    z_chi_denominator = np.sqrt(
        effect_df["chi_sd_Control"] ** 2
        + effect_df["chi_sd_CCl4"] ** 2
    )

    effect_df["Z_chi"] = (
        effect_df["chi_mean_Control"]
        - effect_df["chi_mean_CCl4"]
    ) / z_chi_denominator

    # ========================================================
    # E. Remove invalid effect metrics
    # ========================================================

    effect_df = (
        effect_df
        .replace([np.inf, -np.inf], np.nan)
        .dropna(
            subset=[
                "Z_xi",
                "Z_chi",
                "log_xi_ratio",
                "delta_chi",
            ]
        )
        .copy()
    )

    print(
        f"After valid paper metrics: "
        f"{len(effect_df)} genes"
    )

    # ========================================================
    # F. Final column order
    # ========================================================

    final_columns = [
        "Gene",
        "xi_mean_Control",
        "xi_mean_CCl4",
        "xi_sd_Control",
        "xi_sd_CCl4",
        "xi_n_Control",
        "xi_n_CCl4",
        "r2_mean_Control",
        "r2_mean_CCl4",
        "chi_mean_Control",
        "chi_mean_CCl4",
        "chi_sd_Control",
        "chi_sd_CCl4",
        "chi_n_Control",
        "chi_n_CCl4",
        "Z_xi",
        "Z_chi",
        "log_xi_ratio",
        "delta_chi",
    ]

    effect_df = effect_df[final_columns].copy()

    print("\n==============================================")
    print("Stage 3 Complete")
    print("==============================================")
    print("\nCondition-level effect table:")
    print(effect_df)

    return effect_df


if __name__ == "__main__":
    print(
        "Import summarize_effects() into the pipeline driver "
        "and pass it the Stage 2 results."
    )
