#!/usr/bin/env python3
"""Phase 5 — scVI training and deconvolution cross-check."""
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

warnings.filterwarnings("ignore")
SEED = 42
scvi.settings.seed = SEED

ROOT = Path(__file__).resolve().parent.parent


def run_scvi():
    adata = sc.read_h5ad(ROOT / "objects" / "sc_wt_annotated.h5ad")
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="replicate",
        labels_key="cell_type",
    )
    model = scvi.model.SCVI(adata, n_latent=20, gene_likelihood="nb")
    model.train(max_epochs=50, early_stopping=True, plan_kwargs={"lr": 1e-3})
    model_dir = ROOT / "objects" / "scvi_model"
    model.save(model_dir, overwrite=True)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)

    sc.pl.umap(adata, color=["cell_type", "replicate"], show=False)
    plt.savefig(ROOT / "figures" / "scvi_umap.png", bbox_inches="tight")
    plt.savefig(ROOT / "figures" / "scvi_umap.pdf", bbox_inches="tight")
    plt.close()

    # Cell-type signature proxy vs bulk (Spearman on normalized expression)
    norm = model.get_normalized_expression(library_size=1e4)
    if hasattr(norm, "toarray"):
        norm = norm.toarray()
    ct_df = pd.DataFrame(norm, index=adata.obs_names)
    ct_means = ct_df.groupby(adata.obs["cell_type"].values).mean()

    bulk = pd.read_csv(ROOT / "tables" / "bulk_normalized_counts.csv", index_col=0)
    shared = ct_means.columns.intersection(bulk.index)
    ct_means = ct_means[shared]
    bulk_mean = bulk.loc[shared].mean(axis=1)

    scvi_props = {}
    for ct in ct_means.index:
        rho = pd.Series(ct_means.loc[ct], index=shared).corr(
            pd.Series(bulk_mean, index=shared), method="spearman"
        )
        scvi_props[ct] = max(float(rho), 0.0)
    total = sum(scvi_props.values()) or 1.0
    scvi_prop_df = pd.DataFrame({
        "cell_type": list(scvi_props.keys()),
        "scvi_proxy_proportion": [v / total for v in scvi_props.values()],
    })
    scvi_prop_df.to_csv(ROOT / "tables" / "scvi_proxy_proportions.csv", index=False)

    bisque_path = ROOT / "tables" / "bisque_deconvolution_proportions.csv"
    if bisque_path.exists():
        bisque = pd.read_csv(bisque_path, index_col=0)
        bisque_mean = bisque.mean(axis=1)
        bisque_mean = bisque_mean / bisque_mean.sum()
        compare = pd.DataFrame({
            "cell_type": bisque_mean.index,
            "bisque_mean": bisque_mean.values,
        }).merge(scvi_prop_df, on="cell_type", how="outer").fillna(0)
        compare["abs_diff"] = (compare["bisque_mean"] - compare["scvi_proxy_proportion"]).abs()
        compare.to_csv(ROOT / "tables" / "deconvolution_bisque_vs_scvi.csv", index=False)

    adata.write_h5ad(ROOT / "objects" / "sc_wt_annotated.h5ad")
    print("Phase 5 (scVI) complete.")


if __name__ == "__main__":
    run_scvi()
