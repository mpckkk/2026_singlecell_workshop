#!/usr/bin/env python3
"""Build annotated h5ad from Seurat export + Phase 2 PAGA trajectory."""
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, facecolor="white")

ROOT = Path(__file__).resolve().parent.parent
SEED = 42
np.random.seed(SEED)


def build_h5ad():
    exp_dir = ROOT / "objects" / "h5ad_export"
    if not (exp_dir / "matrix.mtx").exists():
        raise FileNotFoundError("Run Phase 1 (Seurat) first to create h5ad_export/")

    X = scipy.io.mmread(exp_dir / "matrix.mtx").T.tocsr()  # cells x genes
    barcodes = pd.read_csv(exp_dir / "barcodes.tsv", header=None)[0].values
    features = pd.read_csv(exp_dir / "features.tsv", header=None)[0].values
    meta = pd.read_csv(exp_dir / "meta.csv", index_col=0)
    umap = pd.read_csv(exp_dir / "umap.csv", index_col=0)
    pca = pd.read_csv(exp_dir / "pca.csv", index_col=0)

    adata = sc.AnnData(X=X, obs=meta.loc[barcodes], var=pd.DataFrame(index=features))
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["X_umap"] = umap.loc[barcodes].values
    adata.obsm["X_pca"] = pca.loc[barcodes].values
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype(str)
    adata.write_h5ad(ROOT / "objects" / "sc_wt_annotated.h5ad")
    return adata


def run_paga(adata):
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=20, n_pcs=20, random_state=SEED)
    sc.tl.paga(adata, groups="cell_type")
    sc.pl.paga(adata, color="cell_type", show=False)
    plt.savefig(ROOT / "figures" / "paga_graph.png", bbox_inches="tight")
    plt.savefig(ROOT / "figures" / "paga_graph.pdf", bbox_inches="tight")
    plt.close()

    sc.tl.umap(adata, init_pos="paga", random_state=SEED)
    sc.pl.umap(adata, color=["cell_type", "seurat_clusters"], show=False)
    plt.savefig(ROOT / "figures" / "paga_umap_celltype.png", bbox_inches="tight")
    plt.savefig(ROOT / "figures" / "paga_umap_celltype.pdf", bbox_inches="tight")
    plt.close()

    if "QC" in adata.obs["cell_type"].cat.categories:
        sc.tl.diffmap(adata, random_state=SEED)
        root_idx = np.flatnonzero(adata.obs["cell_type"].values == "QC")
        if len(root_idx) > 0:
            adata.uns["iroot"] = int(root_idx[0])
            sc.tl.dpt(adata)
            sc.pl.umap(adata, color="dpt_pseudotime", show=False)
            plt.savefig(ROOT / "figures" / "paga_dpt_pseudotime.png", bbox_inches="tight")
            plt.savefig(ROOT / "figures" / "paga_dpt_pseudotime.pdf", bbox_inches="tight")
            plt.close()

    conn = adata.uns["paga"]["connectivities"].toarray()
    ct = list(adata.obs["cell_type"].cat.categories)
    pd.DataFrame(conn, index=ct, columns=ct).to_csv(ROOT / "tables" / "paga_connectivity.csv")

    adata.write_h5ad(ROOT / "objects" / "sc_wt_annotated.h5ad")
    return adata


def describe_trajectories(adata):
    lines = ["# Trajectory summary (PAGA + DPT)\n"]
    if "dpt_pseudotime" in adata.obs.columns:
        for ct in adata.obs["cell_type"].cat.categories:
            sub = adata.obs.loc[adata.obs["cell_type"] == ct, "dpt_pseudotime"]
            lines.append(f"- **{ct}**: median pseudotime = {sub.median():.3f} (n={len(sub)})")
    else:
        lines.append("- DPT not computed (QC cells absent); PAGA connectivity used for topology.")
    (ROOT / "tables" / "trajectory_summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    adata = build_h5ad()
    adata = run_paga(adata)
    describe_trajectories(adata)
    print("Phase 2 complete.")
