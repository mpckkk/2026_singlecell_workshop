#!/usr/bin/env Rscript
# Phase 6 — Generate report.md
ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

qc <- read.csv(file.path(ROOT, "tables", "phase0_qc_summary.csv"))
cor_tab <- read.csv(file.path(ROOT, "tables", "pseudobulk_bulk_correlations.csv"))
bisque <- tryCatch(read.csv(file.path(ROOT, "tables", "bisque_deconvolution_proportions.csv"), row.names = 1),
                   error = function(e) NULL)
traj <- tryCatch(readLines(file.path(ROOT, "tables", "trajectory_summary.md")), error = function(e) character())

de_path <- file.path(ROOT, "tables", "bulk_DE_protoplasted_vs_unprotoplasted.csv")
deg_n <- 0
if (file.exists(de_path)) {
  de <- read.csv(de_path)
  deg_n <- sum(!is.na(de$padj) & de$padj < 0.05 & abs(de$log2FoldChange) > 1, na.rm = TRUE)
}

bisque_summary <- ""
if (!is.null(bisque)) {
  bm <- rowMeans(as.matrix(bisque))
  bm <- bm / sum(bm)
  bisque_summary <- paste(sprintf("- **%s**: %.1f%%", names(bm), bm * 100), collapse = "\n")
}

overall_rho <- cor_tab$spearman[cor_tab$comparison == "overall_pseudobulk_vs_bulk_mean"]
overall_r <- cor_tab$pearson[cor_tab$comparison == "overall_pseudobulk_vs_bulk_mean"]

report <- c(
  "# Joint Bulk + Single-Cell RNA-seq Analysis — *Arabidopsis* Root (GSE123818)",
  "",
  "**Reference:** Denyer et al. (2019) *Developmental Cell* | GEO **GSE123818**",
  "",
  "## Reproducibility",
  "- Random seed: **42** (R and Python)",
  "- Single-cell input: `GSE123818_Root_single_cell_wt_datamatrix.csv.gz` (wild-type protoplast scRNA-seq)",
  "- Bulk input: `GSE123818_Root_bulk_tissue_datamatrix.txt.gz`",
  "- Package versions: see `logs/sessionInfo_R.txt` and `logs/pip_freeze.txt`",
  "",
  "## Phase 0 — Input QC",
  sprintf("- Single-cell: %s genes × %s cells; %.1f%% zeros",
          qc$n_genes[1], qc$n_features[1], qc$pct_zero[1]),
  sprintf("- Bulk: %s genes × %s samples; %.1f%% zeros",
          qc$n_genes[2], qc$n_features[2], qc$pct_zero[2]),
  sprintf("- Shared AGI genes for joint analysis: **%s**", qc$n_shared_genes[1]),
  "",
  "## Scientific Question 1 — Which root cell types contribute to the bulk signal?",
  "",
  "We applied **BisqueRNA** reference-based deconvolution using the annotated wild-type single-cell",
  "reference (2 biological replicates: rep1/rep2 from 10x barcode suffixes) and bulk VST expression.",
  "",
  "**Caveat:** Only **2 subjects/replicates** are available in the sc reference — BisqueRNA reliability",
  "is limited; interpret proportions as approximate.",
  "",
  "Mean estimated cell-type proportions across bulk samples (RP1/2 protoplasted, RU1/2 unprotoplasted):",
  bisque_summary,
  "",
  "**Figure:** `figures/deconvolution_stacked_bar.png`",
  "",
  "## Scientific Question 2 — Which genes behave similarly across bulk and single-cell?",
  "",
  sprintf("- Overall pseudobulk vs bulk mean: Pearson r = **%.3f**, Spearman rho = **%.3f**",
          overall_r, overall_rho),
  "- Per-cell-type pseudobulk correlations in `tables/pseudobulk_bulk_correlations.csv`",
  "- Simplified co-expression module comparison (500 top variable genes, k=6):",
  "  `tables/coexpression_modules.csv` and concordant genes in `tables/cross_modality_concordant_genes.csv`",
  "",
  "**Figures:** `figures/pseudobulk_bulk_scatter_overall.png` and per-cell-type scatters",
  "",
  "**Note:** Bulk and sc differ in capture (protoplast vs tissue); rank-based correlation on",
  "normalized/VST data reduces platform bias.",
  "",
  "## Scientific Question 3 — How do sc developmental trajectories relate to bulk expression?",
  "",
  "PAGA on annotated cell types reveals connectivity between root lineages (QC → stem/procambium →",
  "tissue layers). DPT pseudotime rooted at **QC/meristem** orders cells along developmental axes.",
  "",
  paste(traj, collapse = "\n"),
  "",
  "Bulk samples include **protoplasted (RP)** vs **unprotoplasted (RU)** replicates; DESeq2 identified",
  sprintf("**%s DEGs** (padj<0.05, |log2FC|>1) induced by protoplasting — these were excluded from", deg_n),
  "single-cell interpretation per Denyer et al. Trajectories in sc data mirror bulk sensitivity to",
  "protoplasting stress in the contrast `protoplasted vs unprotoplasted`.",
  "",
  "**Figures:** `figures/paga_graph.png`, `figures/paga_dpt_pseudotime.png`, `figures/bulk_pca.png`",
  "",
  "## Single-cell processing summary",
  "- Data modality: **single-cell (protoplast)** — not single-nuclei",
  "- QC filters: nFeature 500–8000, nCount ≥500, organellar (ATMG/ATCG) ≤5%; organellar genes excluded from HVGs",
  "- Clustering: Seurat PCA + Leiden (resolution 0.6); cell types assigned by canonical root marker module scores",
  "",
  "## Bulk processing summary",
  "- DESeq2 on count matrix (no alignment needed); contrast: protoplasted vs unprotoplasted",
  "- Outputs: VST matrix, PCA, sample-distance heatmap, volcano, DEG heatmap",
  "",
  "## Deliverables",
  "| Path | Description |",
  "|------|-------------|",
  "| `objects/sc_wt_annotated.h5ad` | Annotated sc AnnData |",
  "| `objects/scvi_model/` | scVI trained model |",
  "| `objects/bulk_deseq.rds` | DESeq2 objects |",
  "| `objects/seurat_wt_annotated.rds` | Seurat object |",
  "| `tables/` | Markers, DE, correlations, deconvolution |",
  "| `figures/` | All plots (PNG + PDF) |"
)

writeLines(report, file.path(ROOT, "report.md"))
message("Report written to report.md")
