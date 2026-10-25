
# ============================================================
# XENIUM GENE–GENE CROSS-CORRELATION PIPELINE
# STAGE 3 — Condition-level Effect Table
#
# Healthy vs Diseased
# ============================================================

import numpy as np
import pandas as pd


# ============================================================
# Parameters
# ============================================================
min_valid_cores_per_condition = 3

# ============================================================
# Main Stage 3 function
# ============================================================

def summarize_effects(xi_batch_pair_df):

    # ========================================================
    # STAGE 3 QC
    # ========================================================

    print("\n==============================================")
    print("STAGE 3 — Condition-level pair statistics")
    print("==============================================")

    initial_unique_pairs = (
        xi_batch_pair_df[
            ["Gene1", "Gene2"]
        ]
        .drop_duplicates()
        .shape[0]
    )

    print(
        f"Initial unique gene pairs: "
        f"{initial_unique_pairs:,}"
    )

    # ========================================================
    # A. ξ summary from accepted per-core ξ fits
    # ========================================================

    xi_valid_df = (
        xi_batch_pair_df
        .dropna(subset=["xi"])
        .copy()
    )

    pairs_with_any_valid_xi = (
        xi_valid_df[
            ["Gene1", "Gene2"]
        ]
        .drop_duplicates()
        .shape[0]
    )

    print(
        f"Pairs with at least one valid ξ fit: "
        f"{pairs_with_any_valid_xi:,}"
    )

    xi_summary = (
        xi_valid_df
        .groupby(
            ["Gene1", "Gene2", "Condition"],
            observed=True
        )
        .agg(
            xi_mean=("xi", "mean"),
            xi_sd=("xi", "std"),
            xi_n=("xi", "count"),
            r2_mean=("r2", "mean")
        )
        .reset_index()
    )

    xi_pivot = xi_summary.pivot(
        index=["Gene1", "Gene2"],
        columns="Condition",
        values=[
            "xi_mean",
            "xi_sd",
            "xi_n",
            "r2_mean"
        ]
    )

    xi_pivot.columns = [
        f"{metric}_{condition}"
        for metric, condition in xi_pivot.columns
    ]

    xi_pivot = xi_pivot.reset_index()

    # ========================================================
    # B. χ summary independent of ξ-fit success
    # ========================================================

    chi_valid_df = (
        xi_batch_pair_df
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["chi"])
        .copy()
    )

    pairs_with_valid_chi = (
        chi_valid_df[
            ["Gene1", "Gene2"]
        ]
        .drop_duplicates()
        .shape[0]
    )

    print(
        f"Pairs with valid χ values: "
        f"{pairs_with_valid_chi:,}"
    )

    chi_summary = (
        chi_valid_df
        .groupby(
            ["Gene1", "Gene2", "Condition"],
            observed=True
        )
        .agg(
            chi_mean=("chi", "mean"),
            chi_sd=("chi", "std"),
            chi_n=("chi", "count")
        )
        .reset_index()
    )

    chi_pivot = chi_summary.pivot(
        index=["Gene1", "Gene2"],
        columns="Condition",
        values=[
            "chi_mean",
            "chi_sd",
            "chi_n"
        ]
    )

    chi_pivot.columns = [
        f"{metric}_{condition}"
        for metric, condition in chi_pivot.columns
    ]

    chi_pivot = chi_pivot.reset_index()

    # ========================================================
    # Absolute χ summary (independent of sign)
    # ========================================================

    chi_abs_valid_df = chi_valid_df.copy()

    chi_abs_valid_df["abs_chi"] = (
        chi_abs_valid_df["chi"].abs()
    )

    chi_abs_summary = (
        chi_abs_valid_df
        .groupby(
            ["Gene1", "Gene2", "Condition"],
            observed=True
        )
        .agg(
            abs_chi_mean=("abs_chi", "mean"),
            abs_chi_sd=("abs_chi", "std"),
            abs_chi_n=("abs_chi", "count")
        )
        .reset_index()
    )

    chi_abs_pivot = chi_abs_summary.pivot(
        index=["Gene1", "Gene2"],
        columns="Condition",
        values=[
            "abs_chi_mean",
            "abs_chi_sd",
            "abs_chi_n"
        ]
    )

    chi_abs_pivot.columns = [
        f"{metric}_{condition}"
        for metric, condition in chi_abs_pivot.columns
    ]

    chi_abs_pivot = chi_abs_pivot.reset_index()

    # ========================================================
    # C. Merge ξ and χ summaries
    # ========================================================

    effect_pair_df = (
        xi_pivot
        .merge(
            chi_pivot,
            on=["Gene1", "Gene2"],
            how="inner"
        )
        .merge(
            chi_abs_pivot,
            on=["Gene1", "Gene2"],
            how="inner"
        )
    )

    # Chi-based pivot alternative from the original pipeline:
    #
    # effect_pair_df = (
    #     chi_pivot
    #     .merge(
    #         chi_abs_pivot,
    #         on=["Gene1", "Gene2"],
    #         how="inner"
    #     )
    #     .merge(
    #         xi_pivot,
    #         on=["Gene1", "Gene2"],
    #         how="left"
    #     )
    # )

    print(
        f"After ξ–χ summary merge: "
        f"{len(effect_pair_df):,} pairs"
    )

    # ========================================================
    # Require required ξ and χ summaries for both conditions
    # ========================================================

    required_columns = [
        "xi_mean_Healthy",
        "xi_mean_Diseased",
        "xi_sd_Healthy",
        "xi_sd_Diseased",
        "xi_n_Healthy",
        "xi_n_Diseased",
        "chi_mean_Healthy",
        "chi_mean_Diseased"
    ]

    # Original alternative requiring only χ summaries:
    #
    # required_columns = [
    #     "chi_mean_Healthy",
    #     "chi_mean_Diseased",
    #     "chi_sd_Healthy",
    #     "chi_sd_Diseased",
    #     "chi_n_Healthy",
    #     "chi_n_Diseased"
    # ]

    effect_pair_df = (
        effect_pair_df
        .dropna(subset=required_columns)
        .copy()
    )

    print(
        f"After requiring Healthy and Diseased summaries: "
        f"{len(effect_pair_df):,} pairs"
    )

    # ========================================================
    # Require positive condition-level mean ξ
    # ========================================================

    effect_pair_df = effect_pair_df[
        (effect_pair_df["xi_mean_Healthy"] > 0)
        &
        (effect_pair_df["xi_mean_Diseased"] > 0)
    ].copy()

    print(
        f"After positive mean ξ filter: "
        f"{len(effect_pair_df):,} pairs"
    )

    # ========================================================
    # Require enough accepted ξ fits in each condition
    # ========================================================

    effect_pair_df = effect_pair_df[
        (
            effect_pair_df["xi_n_Healthy"]
            >= min_valid_cores_per_condition
        )
        &
        (
            effect_pair_df["xi_n_Diseased"]
            >= min_valid_cores_per_condition
        )
    ].copy()

    print(
        f"After requiring ≥{min_valid_cores_per_condition} "
        f"valid cores per condition: "
        f"{len(effect_pair_df):,} pairs"
    )

    # ========================================================
    # Paper metrics
    # ========================================================

    # --------------------------------------------------------
    # Condition-level standardized ξ difference
    # --------------------------------------------------------

    z_denominator = np.sqrt(
        effect_pair_df["xi_sd_Healthy"] ** 2
        +
        effect_pair_df["xi_sd_Diseased"] ** 2
    )

    effect_pair_df["Z_xi"] = (
        effect_pair_df["xi_mean_Healthy"]
        -
        effect_pair_df["xi_mean_Diseased"]
    ) / z_denominator

    # --------------------------------------------------------
    # Scale-independent ξ effect magnitude
    # --------------------------------------------------------

    effect_pair_df["log_xi_ratio"] = np.log(
        effect_pair_df["xi_mean_Healthy"]
        /
        effect_pair_df["xi_mean_Diseased"]
    )

    # --------------------------------------------------------
    # Model-free susceptibility difference
    # --------------------------------------------------------

    effect_pair_df["delta_chi"] = (
        effect_pair_df["chi_mean_Healthy"]
        -
        effect_pair_df["chi_mean_Diseased"]
    )

    # --------------------------------------------------------
    # Magnitude-only susceptibility difference
    # --------------------------------------------------------

    effect_pair_df["delta_abs_chi"] = (
        effect_pair_df["abs_chi_mean_Healthy"]
        -
        effect_pair_df["abs_chi_mean_Diseased"]
    )

    # --------------------------------------------------------
    # Condition-level standardized χ difference
    # --------------------------------------------------------

    z_chi_denominator = np.sqrt(
        effect_pair_df["chi_sd_Healthy"] ** 2
        +
        effect_pair_df["chi_sd_Diseased"] ** 2
    )

    effect_pair_df["Z_chi"] = (
        effect_pair_df["chi_mean_Healthy"]
        -
        effect_pair_df["chi_mean_Diseased"]
    ) / z_chi_denominator

    # --------------------------------------------------------
    # Condition-level standardized |χ| difference
    # --------------------------------------------------------

    z_abs_chi_denominator = np.sqrt(
        effect_pair_df["abs_chi_sd_Healthy"] ** 2
        +
        effect_pair_df["abs_chi_sd_Diseased"] ** 2
    )

    effect_pair_df["Z_abs_chi"] = (
        effect_pair_df["abs_chi_mean_Healthy"]
        -
        effect_pair_df["abs_chi_mean_Diseased"]
    ) / z_abs_chi_denominator

    # ========================================================
    # Remove non-finite effect metrics
    # ========================================================

    effect_pair_df = (
        effect_pair_df
        .replace([np.inf, -np.inf], np.nan)
        .dropna(
            subset=[
                "Z_xi",
                "Z_chi",
                "Z_abs_chi",
                "log_xi_ratio",
                "delta_chi",
                "delta_abs_chi"
            ]
        )
        .copy()
    )

    print(
        f"After retaining finite paper metrics: "
        f"{len(effect_pair_df):,} pairs"
    )

    # ========================================================
    # Pair label
    # ========================================================

    effect_pair_df["Pair"] = (
        effect_pair_df["Gene1"].astype(str)
        +
        "–"
        +
        effect_pair_df["Gene2"].astype(str)
    )

    # ========================================================
    # Final column order
    # ========================================================

    effect_pair_df = effect_pair_df[
        [
            "Gene1",
            "Gene2",
            "Pair",

            "xi_mean_Healthy",
            "xi_mean_Diseased",

            "xi_sd_Healthy",
            "xi_sd_Diseased",

            "xi_n_Healthy",
            "xi_n_Diseased",

            "r2_mean_Healthy",
            "r2_mean_Diseased",

            "chi_mean_Healthy",
            "chi_mean_Diseased",

            "chi_sd_Healthy",
            "chi_sd_Diseased",

            "chi_n_Healthy",
            "chi_n_Diseased",

            "abs_chi_mean_Healthy",
            "abs_chi_mean_Diseased",

            "abs_chi_sd_Healthy",
            "abs_chi_sd_Diseased",

            "Z_xi",
            "Z_chi",
            "Z_abs_chi",

            "log_xi_ratio",

            "delta_chi",
            "delta_abs_chi"
        ]
    ].copy()


    return effect_pair_df
