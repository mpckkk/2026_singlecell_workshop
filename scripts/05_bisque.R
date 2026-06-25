#!/usr/bin/env Rscript
# Phase 5 — BisqueRNA deconvolution
set.seed(42)
ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(BisqueRNA)
  library(Biobase)
  library(Seurat)
  library(ggplot2)
  library(SummarizedExperiment)
})

inp <- readRDS(file.path(ROOT, "objects", "phase0_input.rds"))
seu <- readRDS(file.path(ROOT, "objects", "seurat_wt_annotated.rds"))
bulk_obj <- readRDS(file.path(ROOT, "objects", "bulk_deseq.rds"))

shared <- Reduce(intersect, list(inp$shared_genes, rownames(bulk_obj$vsd), rownames(seu)))
counts <- GetAssayData(seu, layer = "counts")
meta <- seu@meta.data

# Reference ExpressionSet (log2 CPM normalized per cell, BisqueRNA convention)
expr_mat <- as.matrix(counts[shared, ])
expr_norm <- apply(expr_mat, 2, function(x) log2(x / sum(x) * 1e6 + 1))

pheno <- data.frame(
  cellType = meta$cell_type,
  SubjectName = meta$SubjectName,
  row.names = rownames(meta)
)
ref_eset <- ExpressionSet(assayData = assayDataNew(exprs = expr_norm),
                           phenoData = AnnotatedDataFrame(pheno))

# Bulk ExpressionSet — map bulk replicates to sc SubjectName (rep1/rep2)
bulk_vst <- SummarizedExperiment::assay(bulk_obj$vsd)[shared, , drop = FALSE]
subject_map <- c(RP1 = "rep1", RP2 = "rep2", RU1 = "rep1", RU2 = "rep2")
bulk_pheno <- data.frame(
  sample = colnames(bulk_vst),
  SubjectName = unname(subject_map[colnames(bulk_vst)]),
  row.names = colnames(bulk_vst)
)
bulk_eset <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(bulk_vst)),
                           phenoData = AnnotatedDataFrame(bulk_pheno))

n_subjects <- length(unique(meta$SubjectName))
message("BisqueRNA reference: ", length(unique(meta$cell_type)), " cell types, ",
        n_subjects, " subjects/replicates")

# Run decomposition
prop <- ReferenceBasedDecomposition(bulk.eset = bulk_eset, sc.eset = ref_eset,
                                    markers = NULL, use.overlap = FALSE)

prop_mat <- prop$bulk.props
write.csv(prop_mat, file.path(ROOT, "tables", "bisque_deconvolution_proportions.csv"))

# Stacked bar plot
prop_long <- as.data.frame(prop_mat)
prop_long$cell_type <- rownames(prop_long)
id_vars <- setdiff(colnames(prop_long), "cell_type")
prop_long <- reshape(
  prop_long,
  direction = "long",
  varying = list(id_vars),
  v.names = "proportion",
  timevar = "sample",
  idvar = "cell_type"
)

p_stack <- ggplot(prop_long, aes(sample, proportion, fill = cell_type)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Bulk sample", y = "Estimated cell-type proportion", fill = "Cell type",
       title = "BisqueRNA deconvolution of bulk root samples") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(ROOT, "figures", "deconvolution_stacked_bar.png"), p_stack, width = 9, height = 6, dpi = 150)
ggsave(file.path(ROOT, "figures", "deconvolution_stacked_bar.pdf"), p_stack, width = 9, height = 6)

saveRDS(list(proportions = prop_mat, n_subjects = n_subjects, ref_pheno = pheno),
        file.path(ROOT, "objects", "bisque_deconvolution.rds"))
message("Phase 5 (BisqueRNA) complete.")
