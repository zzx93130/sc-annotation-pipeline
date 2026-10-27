"""Helpers for sc-annotation-pipeline: the mandatory regress_out -> scale -> PCA path and its checks."""
import numpy as np


def regress_scale_pca(adata, keys=None, n_pcs=50, max_value=10.0, block_genes=500,
                      layer=None, seed=0, n_jobs=8, key_added="X_pca"):
    """Default PCA path for ALL scanpy clustering: sc.pp.regress_out(keys) -> sc.pp.scale(max_value) -> sc.pp.pca(arpack).

    adata : cells x HVG AnnData; X (or `layer`) log-normalised. Not modified except obsm[key_added] and
            uns['regress_scale_pca'].
    keys  : obs columns to regress; default ['total_counts', 'pct_counts_mt'] computed by
            sc.pp.calculate_qc_metrics on raw counts of ALL genes (before HVG subsetting). Pass [] to disable
            (only for a deliberate unregressed control).
    regress_out fits each gene independently, so blocks of `block_genes` give the same result as one call
    with bounded memory. Returns adata.
    """
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp
    if keys is None:
        keys = ["total_counts", "pct_counts_mt"]
    keys = list(keys)
    miss = [k for k in keys if k not in adata.obs.columns]
    if miss:
        raise KeyError(f"regress keys missing from obs: {miss}; run sc.pp.calculate_qc_metrics on raw counts first")
    X = adata.layers[layer] if layer is not None else adata.X
    n, g = X.shape
    cov = pd.DataFrame(adata.obs[keys].to_numpy(dtype=np.float64), columns=keys, index=np.arange(n).astype(str))
    if keys and not np.isfinite(cov.to_numpy()).all():
        raise ValueError("regress covariates contain NaN/inf")
    Xc = X.tocsc() if sp.issparse(X) else np.asarray(X)
    Y = np.empty((n, g), dtype=np.float32)
    for j0 in range(0, g, block_genes):
        j1 = min(g, j0 + block_genes)
        blk = Xc[:, j0:j1]
        blk = blk.toarray() if sp.issparse(blk) else np.array(blk)
        if keys:
            b = ad.AnnData(X=blk.astype(np.float32), obs=cov)
            sc.pp.regress_out(b, keys, n_jobs=n_jobs)
            blk = b.X
        Y[:, j0:j1] = blk
    t = ad.AnnData(X=Y)
    sc.pp.scale(t, max_value=max_value)
    nc = min(n_pcs, min(t.shape) - 1)
    sc.pp.pca(t, n_comps=nc, svd_solver="arpack", random_state=seed, zero_center=True)
    adata.obsm[key_added] = np.asarray(t.obsm["X_pca"], dtype=np.float32)
    adata.uns["regress_scale_pca"] = {"keys": keys, "n_pcs": int(nc), "max_value": float(max_value),
                                      "block_genes": int(block_genes), "seed": int(seed),
                                      "variance_ratio": np.asarray(t.uns["pca"]["variance_ratio"]).tolist()}
    return adata


def qc_pc_correlation(adata, keys=None, rep="X_pca", n_top=10):
    """Max |Pearson r| between the first n_top components of obsm[rep] and each QC covariate
    (total_counts is log1p-transformed). Use before/after regression as the check that it worked."""
    import pandas as pd
    if keys is None:
        keys = ["total_counts", "pct_counts_mt"]
    Z = np.asarray(adata.obsm[rep])[:, :n_top]
    out = {}
    for k in keys:
        v = adata.obs[k].to_numpy(dtype=float)
        if k == "total_counts":
            v = np.log1p(v)
        r = [abs(np.corrcoef(Z[:, i], v)[0, 1]) for i in range(Z.shape[1])]
        out[k] = float(np.nanmax(r))
    return pd.Series(out, name=f"max_abs_r_top{n_top}")


def compartment_drift(obs, before_col, after_col, depth_col="total_counts", batch_col=None):
    """Cells whose compartment/lineage call changes between an unregressed and a regressed run
    (or vs a manual reference). Regressing total_counts can push very shallow cells across lineages;
    report the count, their median depth vs overall, and the batches they come from."""
    import pandas as pd
    a = obs[before_col].astype(str)
    b = obs[after_col].astype(str)
    m = obs[a != b]
    res = {"n_changed": int(len(m)), "frac_changed": float(len(m) / max(len(obs), 1)),
           "median_depth_changed": float(m[depth_col].median()) if len(m) else float("nan"),
           "median_depth_all": float(obs[depth_col].median())}
    if batch_col is not None and len(m):
        res["by_batch"] = m[batch_col].astype(str).value_counts().head(10).to_dict()
    res["transitions"] = (pd.crosstab(a[a != b], b[a != b]).stack().sort_values(ascending=False)
                          .head(10).to_dict() if len(m) else {})
    return res
