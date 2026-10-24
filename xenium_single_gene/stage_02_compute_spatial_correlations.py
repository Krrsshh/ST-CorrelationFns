
# ============================================================
# PIPELINE 5 (XENIUM)
# STAGE 2 — Gene-wise ξ and χ for Every Core
# ============================================================

import numpy as np
import pandas as pd

from scipy import sparse
from scipy.optimize import curve_fit
from sklearn.neighbors import KDTree


# ============================================================
# Spatial parameters
# ============================================================

MAX_DIST = 100          # µm
BIN_SIZE = 5            # µm
FIT_START_R = 30        # µm

MIN_FIT_POINTS = 5
MIN_ABS_CORR = 0.0      # 0.0 disables this filter; use 0.02 to enable it
MIN_HEPATOCYTES = 100
FLAT_CURVE_STD = 1e-4
MIN_XI_R2 = 0.1
MAX_XI_MULTIPLIER = 5


# ============================================================
# Helper — convert sparse expression matrices to dense arrays
# ============================================================

def to_dense(X):
    """Convert a sparse or dense expression matrix to a NumPy array."""
    if sparse.issparse(X):
        return X.toarray()
    return np.asarray(X)


# ============================================================
# Exponential model
# ============================================================

def exp_decay(r, A, xi):
    """Single-exponential decay model: C(r) = A * exp(-r / xi)."""
    return A * np.exp(-r / xi)


# ============================================================
# Helper — build unique spatial pairs once per core
# ============================================================

def build_spatial_pairs(
    coords,
    bins,
    max_dist=MAX_DIST,
):
    """
    Identify unique unordered cell pairs within max_dist.

    Each pair is assigned to a distance bin. Pairs are generated
    once per core and reused across all genes in that core.
    """
    n_bins = len(bins) - 1
    tree = KDTree(coords)

    pair_i_list = []
    pair_j_list = []
    pair_bin_list = []

    for i, coord in enumerate(coords):
        inds = tree.query_radius([coord], r=max_dist)[0]

        # Keep each unordered pair only once.
        inds = inds[inds > i]

        if len(inds) == 0:
            continue

        dists = np.linalg.norm(coords[inds] - coord, axis=1)
        pair_bins = np.searchsorted(bins, dists) - 1

        valid = (pair_bins >= 0) & (pair_bins < n_bins)

        inds = inds[valid]
        pair_bins = pair_bins[valid]

        if len(inds) == 0:
            continue

        pair_i_list.append(
            np.full(len(inds), i, dtype=np.int32)
        )
        pair_j_list.append(inds.astype(np.int32))
        pair_bin_list.append(pair_bins.astype(np.int16))

    if not pair_i_list:
        return (
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int16),
        )

    return (
        np.concatenate(pair_i_list),
        np.concatenate(pair_j_list),
        np.concatenate(pair_bin_list),
    )


# ============================================================
# Helper — compute one gene correlation curve
# ============================================================

def compute_gene_curve(
    x,
    pair_i,
    pair_j,
    pair_bins,
    n_bins,
):
    """
    Compute mean covariance and normalized correlation per
    distance bin for one gene within one core.
    """
    x = np.asarray(x, dtype=np.float64)

    x_centered = x - np.mean(x)
    gene_var = np.var(x)

    if not np.isfinite(gene_var) or gene_var <= 0:
        return None

    pair_products = x_centered[pair_i] * x_centered[pair_j]

    # Sum centered expression products per distance bin.
    bin_sums = np.bincount(
        pair_bins,
        weights=pair_products,
        minlength=n_bins,
    ).astype(float)

    # Count cell pairs per distance bin.
    pair_counts = np.bincount(
        pair_bins,
        minlength=n_bins,
    ).astype(int)

    mean_cov = np.full(n_bins, np.nan, dtype=float)

    nonempty = pair_counts > 0
    mean_cov[nonempty] = bin_sums[nonempty] / pair_counts[nonempty]

    corr_vals = mean_cov / gene_var

    return mean_cov, corr_vals, pair_counts, gene_var


# ============================================================
# Main computation
# ============================================================

def compute_spatial_correlations(
    core_data,
    hep_hvg_genes,
    expression_layer="lognorm",
    max_dist=MAX_DIST,
    bin_size=BIN_SIZE,
    fit_start_r=FIT_START_R,
    min_fit_points=MIN_FIT_POINTS,
    min_abs_corr=MIN_ABS_CORR,
    min_hepatocytes=MIN_HEPATOCYTES,
    flat_curve_std=FLAT_CURVE_STD,
    min_xi_r2=MIN_XI_R2,
    max_xi_multiplier=MAX_XI_MULTIPLIER,
):
    """
    Compute core-level gene-wise spatial correlation statistics.

    Parameters
    ----------
    core_data : dict
        Output from Stage 1. Each core entry contains condition,
        sample_id, and an AnnData object under the key 'adata'.

    hep_hvg_genes : list
        Hepatocyte-specific highly variable genes selected in Stage 1.

    expression_layer : str
        AnnData layer containing the expression values used for
        correlation calculations.

    Returns
    -------
    xi_core_gene_df : pandas.DataFrame
        Core-level gene-wise ξ, χ, fit statistics, and metadata.

    qc : dict
        Counters describing processing and filtering outcomes.
    """

    # --------------------------------------------------------
    # Set up distance bins and fitting/integration region
    # --------------------------------------------------------

    bins = np.arange(0, max_dist + bin_size, bin_size)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins = len(bin_centers)

    # Use the same distance region for ξ fitting and χ integration.
    tail_mask = bin_centers >= fit_start_r

    # --------------------------------------------------------
    # QC counters
    # --------------------------------------------------------

    qc = {
        "initial_genes": len(hep_hvg_genes),
        "processed_gene_core_rows": 0,
        "zero_variance": 0,
        "no_finite_corr": 0,
        "flat_curve": 0,
        "too_few_tail_points": 0,
        "negative_tail": 0,  # Positive-tail filter is disabled.
        "fit_failure": 0,
        "failed_stability": 0,
        "valid_xi": 0,
    }

    all_results = []

    print("\n==============================================")
    print("Computing gene-wise ξ and χ for every core")
    print("==============================================")

    # --------------------------------------------------------
    # Process each core
    # --------------------------------------------------------

    for core_id, info in core_data.items():

        condition = info["condition"]
        sample_id = info["sample_id"]
        sub = info["adata"]

        print("\n----------------------------------------------")
        print(f"Processing {condition}: {core_id} | {sample_id}")
        print("----------------------------------------------")

        n_cells = sub.n_obs

        if n_cells < min_hepatocytes:
            print(f"Skipping: only {n_cells} hepatocytes")
            continue

        if "spatial" not in sub.obsm:
            raise KeyError(
                f"Core {core_id} has no 'spatial' coordinates in obsm."
            )

        if expression_layer not in sub.layers:
            raise KeyError(
                f"Expression layer '{expression_layer}' is missing "
                f"from core {core_id}."
            )

        coords = np.asarray(sub.obsm["spatial"])

        # Expression matrix: cells × hepatocyte-specific HVGs.
        X = sub[:, hep_hvg_genes].layers[expression_layer]
        X = to_dense(X).astype(np.float32)
        gene_names = np.asarray(hep_hvg_genes)

        print(f"HEP cells: {n_cells}")
        print(f"Genes: {X.shape[1]}")

        # ----------------------------------------------------
        # Build spatial pairs once for this core
        # ----------------------------------------------------

        pair_i, pair_j, pair_bins = build_spatial_pairs(
            coords=coords,
            bins=bins,
            max_dist=max_dist,
        )

        print(
            f"Unique cell pairs within {max_dist} µm: "
            f"{len(pair_i):,}"
        )

        if len(pair_i) == 0:
            print("No valid spatial pairs. Skipping.")
            continue

        # ----------------------------------------------------
        # Loop over genes
        # ----------------------------------------------------

        for gene_idx, gene in enumerate(gene_names):

            qc["processed_gene_core_rows"] += 1
            x = X[:, gene_idx]

            curve_out = compute_gene_curve(
                x=x,
                pair_i=pair_i,
                pair_j=pair_j,
                pair_bins=pair_bins,
                n_bins=n_bins,
            )

            if curve_out is None:
                qc["zero_variance"] += 1
                continue

            mean_cov, corr_vals, pair_counts, gene_var = curve_out

            finite_corr = corr_vals[np.isfinite(corr_vals)]

            if len(finite_corr) == 0:
                qc["no_finite_corr"] += 1
                continue

            # Optional filter for curves with very small correlations.
            # In the current pipeline, min_abs_corr=0.0 disables it.
            if np.nanmax(np.abs(finite_corr)) < min_abs_corr:
                qc["flat_curve"] += 1
                continue

            # ------------------------------------------------
            # Susceptibility χ
            # Integrate over the same tail used for ξ fitting.
            # ------------------------------------------------

            chi = np.nansum(
                corr_vals[tail_mask]
                * (2 * np.pi * bin_centers[tail_mask])
                * bin_size
            )

            # ------------------------------------------------
            # Prepare tail for exponential fitting
            # ------------------------------------------------

            r_fit = bin_centers[tail_mask]
            y_fit = corr_vals[tail_mask]

            valid_fit = np.isfinite(r_fit) & np.isfinite(y_fit)
            r_fit = r_fit[valid_fit]
            y_fit = y_fit[valid_fit]

            xi = np.nan
            xi_err = np.nan
            r2 = np.nan
            A_fit = np.nan

            if len(r_fit) < min_fit_points:
                qc["too_few_tail_points"] += 1
                continue

            # Skip essentially flat correlation curves.
            # χ is not stored for these rows, matching the original
            # pipeline's behavior.
            if np.nanstd(corr_vals) < flat_curve_std:
                qc["flat_curve"] += 1
                continue

            # The positive-tail sanity check remains disabled,
            # matching the original notebook.
            # if y_fit[0] <= 0:
            #     qc["negative_tail"] += 1
            #     continue

            # ------------------------------------------------
            # Fit the exponential decay model
            # ------------------------------------------------

            try:
                popt, pcov = curve_fit(
                    exp_decay,
                    r_fit,
                    y_fit,
                    p0=[y_fit[0], 20],
                    maxfev=10000,
                )

                A_candidate, xi_candidate = popt

                xi_error_candidate = (
                    np.sqrt(pcov[1, 1])
                    if pcov is not None
                    else np.nan
                )

                y_pred = exp_decay(
                    r_fit,
                    A_candidate,
                    xi_candidate,
                )

                ss_res = np.sum((y_fit - y_pred) ** 2)
                ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)

                r2_candidate = (
                    1 - ss_res / ss_tot
                    if ss_tot > 0
                    else 0
                )

                # Accept ξ only when all stability criteria pass.
                if (
                    xi_candidate > 0
                    and xi_candidate < max_dist * max_xi_multiplier
                    and r2_candidate > min_xi_r2
                ):
                    xi = xi_candidate
                    xi_err = xi_error_candidate
                    r2 = r2_candidate
                    A_fit = A_candidate
                    qc["valid_xi"] += 1

                else:
                    # Keep χ and the fitted R², but reject ξ.
                    r2 = r2_candidate
                    qc["failed_stability"] += 1

            except Exception:
                qc["fit_failure"] += 1

            # ------------------------------------------------
            # Store result
            # ------------------------------------------------

            all_results.append({
                "Core": core_id,
                "Sample": sample_id,
                "Condition": condition,
                "Gene": gene,
                "n_hepatocytes": n_cells,
                "n_spatial_pairs": len(pair_i),
                "xi": xi,
                "xi_err": xi_err,
                "chi": chi,
                "A_fit": A_fit,
                "r2": r2,
                "gene_variance": gene_var,
            })

        print(f"Completed genes for {core_id}: {len(gene_names)}")

    # --------------------------------------------------------
    # Assemble per-core results table
    # --------------------------------------------------------

    result_columns = [
        "Core",
        "Sample",
        "Condition",
        "Gene",
        "n_hepatocytes",
        "n_spatial_pairs",
        "xi",
        "xi_err",
        "chi",
        "A_fit",
        "r2",
        "gene_variance",
    ]

    xi_core_gene_df = pd.DataFrame(
        all_results,
        columns=result_columns,
    )

    print("\n==============================================")
    print("Stage 2 Complete")
    print("==============================================")
    print("Result dimensions:", xi_core_gene_df.shape)

    print("\nPer-core gene-wise results:")
    print(xi_core_gene_df.head())

    print("\nValid ξ estimates by condition:")
    if not xi_core_gene_df.empty:
        print(
            xi_core_gene_df.groupby("Condition")["xi"]
            .agg(Valid_xi="count", Total_gene_core_rows="size")
        )
    else:
        print("No gene-core result rows were generated.")

    print("\n========================================")
    print("Stage 2 QC")
    print("========================================")

    for key, value in qc.items():
        print(f"{key:30s}: {value}")

    return xi_core_gene_df, qc
