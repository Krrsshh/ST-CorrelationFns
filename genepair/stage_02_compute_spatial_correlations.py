
# ============================================================
# XENIUM GENE–GENE CROSS-CORRELATION PIPELINE
# STAGE 2 — Spatial Correlations, ξ Fitting, χ and QC
#
# Requires Stage 1 objects:
#   core_data
#   genes_use
# ============================================================

import numpy as np
import pandas as pd

from scipy import sparse
from scipy.optimize import curve_fit
from sklearn.neighbors import KDTree


# ============================================================
# Spatial and fitting parameters
# ============================================================

max_dist = 100
bin_size = 5
fit_start_r = 30

min_fit_points = 5
min_abs_corr = 0.0
min_pairs_per_bin = 1


# ============================================================
# Exponential decay model
# ============================================================

def exp_decay(r, A, xi):

    return A * np.exp(-r / xi)


# ============================================================
# Helper — convert sparse matrices to dense arrays
# ============================================================

def to_dense(X):

    if sparse.issparse(X):
        return X.toarray()

    return np.asarray(X)


# ============================================================
# Helper — build all unique spatial pairs once per core
# ============================================================

def build_spatial_pairs(coords, bins, n_bins):

    tree = KDTree(coords)

    pair_i_list = []
    pair_j_list = []
    pair_bin_list = []

    for i, coord in enumerate(coords):

        inds = tree.query_radius(
            [coord],
            r=max_dist
        )[0]

        inds = inds[inds > i]

        if len(inds) == 0:
            continue

        dists = np.linalg.norm(
            coords[inds] - coord,
            axis=1
        )

        pair_bins = (
            np.searchsorted(bins, dists) - 1
        )

        valid = (
            (pair_bins >= 0)
            &
            (pair_bins < n_bins)
        )

        inds = inds[valid]
        pair_bins = pair_bins[valid]

        if len(inds) == 0:
            continue

        pair_i_list.append(
            np.full(
                len(inds),
                i,
                dtype=np.int32
            )
        )

        pair_j_list.append(
            inds.astype(np.int32)
        )

        pair_bin_list.append(
            pair_bins.astype(np.int16)
        )

    if len(pair_i_list) == 0:

        return (
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int16)
        )

    pair_i = np.concatenate(pair_i_list)
    pair_j = np.concatenate(pair_j_list)
    pair_bins = np.concatenate(pair_bin_list)

    return pair_i, pair_j, pair_bins


# ============================================================
# Helper — compute all gene-pair correlation curves
# ============================================================

def compute_all_pair_curves(
    X,
    pair_i,
    pair_j,
    pair_bins,
    n_genes,
    n_bins
):

    X = np.asarray(
        X,
        dtype=np.float64
    )

    gene_means = X.mean(axis=0)
    X_centered = X - gene_means

    gene_variances = X.var(axis=0)
    gene_sds = np.sqrt(gene_variances)

    denominator = np.outer(
        gene_sds,
        gene_sds
    )

    denominator[
        (~np.isfinite(denominator))
        |
        (denominator <= 0)
    ] = np.nan

    corr_bins = np.full(
        (n_bins, n_genes, n_genes),
        np.nan,
        dtype=np.float64
    )

    spatial_pair_counts = np.bincount(
        pair_bins,
        minlength=n_bins
    ).astype(np.int64)

    for bin_idx in range(n_bins):

        in_bin = (
            pair_bins == bin_idx
        )

        n_pairs_bin = int(
            spatial_pair_counts[bin_idx]
        )

        if n_pairs_bin < min_pairs_per_bin:
            continue

        cell_i = pair_i[in_bin]
        cell_j = pair_j[in_bin]

        Xi = X_centered[cell_i]
        Xj = X_centered[cell_j]

        directed_sum = Xi.T @ Xj

        symmetric_sum = (
            0.5
            *
            (
                directed_sum
                +
                directed_sum.T
            )
        )

        mean_covariance = (
            symmetric_sum
            /
            n_pairs_bin
        )

        corr_bins[bin_idx] = (
            mean_covariance
            /
            denominator
        )

    return (
        corr_bins,
        spatial_pair_counts,
        gene_variances
    )


# ============================================================
# Main Stage 2 function
# ============================================================

def compute_spatial_correlations(core_data, genes_use):

    # ========================================================
    # Distance bins and gene-pair indices
    # ========================================================

    bins = np.arange(
        0,
        max_dist + bin_size,
        bin_size
    )

    bin_centers = (
        bins[:-1] + bins[1:]
    ) / 2

    n_bins = len(bin_centers)

    # Same distance region for ξ fitting and χ integration
    tail_mask = bin_centers >= fit_start_r

    n_genes = len(genes_use)

    gene1_idx, gene2_idx = np.triu_indices(
        n_genes,
        k=1
    )

    n_gene_pairs = len(gene1_idx)

    gene_names_array = np.asarray(
        genes_use,
        dtype=object
    )

    pair_gene1_names = gene_names_array[gene1_idx]
    pair_gene2_names = gene_names_array[gene2_idx]

    print("\n==============================================")
    print("Stage 2 Configuration")
    print("==============================================")

    print(f"Genes: {n_genes}")
    print(f"Gene pairs per core: {n_gene_pairs:,}")
    print(f"Cores: {len(core_data)}")

    print(
        "Total core–gene-pair combinations:",
        f"{len(core_data) * n_gene_pairs:,}"
    )

    # ========================================================
    # Storage initialized for Stage 2B
    # ========================================================

    all_pair_results = []
    core_paircount_results = {}
    core_qc_rows = []

    # ========================================================
    # Stage 2B — Main computation, ξ fitting, χ and QC
    # ========================================================

    print("\n==============================================")
    print("Computing gene-pair ξ and χ for every batch")
    print("==============================================")

    for core_id, info in core_data.items():

        condition = info["condition"]
        sample_id = info["sample_id"]
        adata = info["adata"]

        print("\n----------------------------------------------")
        print(f"Processing {condition}: {core_id}")
        print("----------------------------------------------")

        n_hepatocytes = adata.n_obs

        if n_hepatocytes < 100:
            print(f"Skipping: only {n_hepatocytes} HEP cells")
            continue

        # ====================================================
        # Coordinates and selected-gene expression
        # ====================================================

        coords = adata.obs[
            ["x_centroid", "y_centroid"]
        ].to_numpy(dtype=np.float64)

        X = to_dense(
            adata.X
        ).astype(
            np.float32,
            copy=False
        )

        if not np.isfinite(X).all():
            raise ValueError(
                f"Non-finite expression values found in {core_id}."
            )

        print(f"HEP cells: {n_hepatocytes:,}")
        print(f"Genes: {X.shape[1]}")

        # ====================================================
        # Build unique spatial pairs once for this core
        # ====================================================

        pair_i, pair_j, pair_bins = build_spatial_pairs(
            coords,
            bins,
            n_bins
        )

        n_spatial_pairs = len(pair_i)

        print(
            f"Unique HEP pairs within {max_dist} µm: "
            f"{n_spatial_pairs:,}"
        )

        if n_spatial_pairs == 0:
            print("No valid spatial pairs. Skipping batch.")
            continue

        # ====================================================
        # Compute all gene-pair curves simultaneously
        # ====================================================

        (
            corr_bins,
            spatial_pair_counts,
            gene_variances
        ) = compute_all_pair_curves(
            X=X,
            pair_i=pair_i,
            pair_j=pair_j,
            pair_bins=pair_bins,
            n_genes=n_genes,
            n_bins=n_bins
        )

        core_paircount_results[core_id] = {
            "Condition": condition,
            "bin_centers": bin_centers.copy(),
            "pair_counts": spatial_pair_counts.copy(),
            "n_spatial_pairs": n_spatial_pairs
        }

        # Extract only the upper-triangle gene pairs
        pair_corr_curves = corr_bins[
            :,
            gene1_idx,
            gene2_idx
        ].T

        pair_var1 = gene_variances[gene1_idx]
        pair_var2 = gene_variances[gene2_idx]

        # ====================================================
        # Batch QC counters
        # ====================================================

        qc = {
            "Core": core_id,
            "Condition": condition,
            "n_hepatocytes": n_hepatocytes,
            "n_spatial_pairs": n_spatial_pairs,
            "initial_gene_pairs": n_gene_pairs,
            "invalid_variance": 0,
            "no_finite_correlation": 0,
            "flat_curve": 0,
            "too_few_tail_points": 0,
            "nonpositive_tail": 0,
            "fit_failure": 0,
            "failed_stability": 0,
            "valid_xi": 0
        }

        batch_rows = []

        # ====================================================
        # Fit ξ for every gene pair
        # ====================================================

        for pair_idx in range(n_gene_pairs):

            gene1 = pair_gene1_names[pair_idx]
            gene2 = pair_gene2_names[pair_idx]

            var1 = pair_var1[pair_idx]
            var2 = pair_var2[pair_idx]

            corr_curve = pair_corr_curves[pair_idx]

            xi = np.nan
            xi_err = np.nan
            chi = np.nan
            A_fit = np.nan
            r2 = np.nan

            filter_status = "Unprocessed"

            # ------------------------------------------------
            # Variance validity
            # ------------------------------------------------

            if (
                not np.isfinite(var1)
                or not np.isfinite(var2)
                or var1 <= 0
                or var2 <= 0
            ):

                qc["invalid_variance"] += 1
                filter_status = "Invalid variance"

            else:

                finite_corr = corr_curve[
                    np.isfinite(corr_curve)
                ]

                # --------------------------------------------
                # Correlation validity
                # --------------------------------------------

                if len(finite_corr) == 0:

                    qc["no_finite_correlation"] += 1
                    filter_status = "No finite correlation"

                elif np.nanmax(
                    np.abs(finite_corr)
                ) < min_abs_corr:

                    qc["flat_curve"] += 1
                    filter_status = "Below correlation threshold"

                else:

                    # ========================================
                    # Susceptibility χ
                    # Same tail region used for ξ fitting
                    # ========================================

                    chi = np.nansum(
                        corr_curve[tail_mask]
                        * (2 * np.pi * bin_centers[tail_mask])
                        * bin_size
                    )

                    # ========================================
                    # Tail data
                    # ========================================

                    r_fit = bin_centers[tail_mask]
                    y_fit = corr_curve[tail_mask]

                    valid_fit = (
                        np.isfinite(r_fit)
                        &
                        np.isfinite(y_fit)
                    )

                    r_fit = r_fit[valid_fit]
                    y_fit = y_fit[valid_fit]

                    # ========================================
                    # Fit eligibility
                    # ========================================

                    if len(r_fit) < min_fit_points:

                        qc["too_few_tail_points"] += 1
                        filter_status = "Too few tail points"

                    # elif y_fit[0] <= 0:
                    #
                    #     qc["nonpositive_tail"] += 1
                    #     filter_status = "Tail starts nonpositive"

                    elif np.nanstd(finite_corr) < 1e-4:

                        qc["flat_curve"] += 1
                        filter_status = "Flat correlation curve"

                    else:

                        try:

                            popt, pcov = curve_fit(
                                exp_decay,
                                r_fit,
                                y_fit,
                                p0=[y_fit[0], 20],
                                maxfev=10000
                            )

                            A_candidate, xi_candidate = popt

                            xi_error_candidate = (
                                np.sqrt(pcov[1, 1])
                                if (
                                    pcov is not None
                                    and np.isfinite(pcov[1, 1])
                                    and pcov[1, 1] >= 0
                                )
                                else np.nan
                            )

                            y_pred = exp_decay(
                                r_fit,
                                A_candidate,
                                xi_candidate
                            )

                            ss_res = np.sum(
                                (y_fit - y_pred) ** 2
                            )

                            ss_tot = np.sum(
                                (y_fit - np.mean(y_fit)) ** 2
                            )

                            r2_candidate = (
                                1 - ss_res / ss_tot
                                if ss_tot > 0
                                else 0
                            )

                            # =================================
                            # Stability filters
                            # =================================

                            if (
                                xi_candidate > 0
                                and xi_candidate < max_dist * 5
                                and r2_candidate > 0.1
                            ):

                                xi = xi_candidate
                                xi_err = xi_error_candidate
                                A_fit = A_candidate
                                r2 = r2_candidate

                                qc["valid_xi"] += 1
                                filter_status = "Accepted"

                            else:

                                r2 = r2_candidate

                                qc["failed_stability"] += 1
                                filter_status = "Failed stability"

                        except Exception:

                            qc["fit_failure"] += 1
                            filter_status = "Fit failure"

            # =================================================
            # Store every pair, including rejected fits
            # =================================================

            batch_rows.append({
                "Core": core_id,
                "Sample": sample_id,
                "Condition": condition,
                "Gene1": gene1,
                "Gene2": gene2,
                "Pair": f"{gene1}–{gene2}",
                "n_hepatocytes": n_hepatocytes,
                "n_spatial_pairs": n_spatial_pairs,
                "var_Gene1": var1,
                "var_Gene2": var2,
                "xi": xi,
                "xi_err": xi_err,
                "chi": chi,
                "A_fit": A_fit,
                "r2": r2,
                "filter_status": filter_status
            })

        # ====================================================
        # Save batch output
        # ====================================================

        all_pair_results.extend(batch_rows)
        core_qc_rows.append(qc)

        print(f"Gene pairs processed: {n_gene_pairs:,}")
        print(f"Accepted ξ fits: {qc['valid_xi']:,}")
        print(f"Nonpositive tails: {qc['nonpositive_tail']:,}")

        print(
            "Failed stability filters:",
            f"{qc['failed_stability']:,}"
        )

    # ========================================================
    # Combined results tables
    # ========================================================

    xi_batch_pair_df = pd.DataFrame(all_pair_results)
    pair_qc_df = pd.DataFrame(core_qc_rows)

    # ========================================================
    # Stage 2 summary
    # ========================================================

    print("\n==============================================")
    print("STAGE 2 COMPLETE")
    print("==============================================")

    print("Result dimensions:", xi_batch_pair_df.shape)

    expected_rows = (
        len(core_qc_rows) * n_gene_pairs
    )

    print("Expected core–pair rows:", expected_rows)
    print("Observed core–pair rows:", len(xi_batch_pair_df))
    print("\nPer-core gene-pair results:")

    print(
        xi_batch_pair_df.head().to_string(index=False)
    )

    print("\nCore-level QC summary:")

    print(
        pair_qc_df.to_string(index=False)
    )

    # ========================================================
    # Overall QC totals
    # ========================================================

    qc_count_columns = [
        "initial_gene_pairs",
        "invalid_variance",
        "no_finite_correlation",
        "flat_curve",
        "too_few_tail_points",
        "nonpositive_tail",
        "fit_failure",
        "failed_stability",
        "valid_xi"
    ]

    overall_qc = pair_qc_df[
        qc_count_columns
    ].sum()

    print("\n==============================================")
    print("Overall Stage 2 QC")
    print("==============================================")

    for metric, value in overall_qc.items():

        print(
            f"{metric:30s}: "
            f"{int(value):,}"
        )

    # ========================================================
    # Valid ξ estimates by condition
    # ========================================================

    print("\nValid ξ estimates by condition:")

    valid_xi_summary = (
        xi_batch_pair_df
        .groupby("Condition")["xi"]
        .agg(["count", "size"])
        .rename(
            columns={
                "count": "Valid_xi",
                "size": "Total_pair_batch_rows"
            }
        )
    )

    print(
        valid_xi_summary.to_string()
    )

    # ========================================================
    # Filter-status QC
    # ========================================================

    print("\nFilter-status counts:")

    filter_status_summary = (
        xi_batch_pair_df["filter_status"]
        .value_counts()
        .rename_axis("filter_status")
        .reset_index(name="count")
    )

    print(
        filter_status_summary.to_string(index=False)
    )

    # ========================================================
    # Return objects required by Stage 3
    # ========================================================

    return (
        xi_batch_pair_df,
        pair_qc_df,
        core_paircount_results
    )
