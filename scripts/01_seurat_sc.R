#!/usr/bin/env Rscript
# Phase 1 — Single-cell processing (Seurat)
set.seed(42)
ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
})

inp <- readRDS(file.path(ROOT, "objects", "phase0_input.rds"))
counts <- inp$sc

# Canonical Arabidopsis root markers (AGI IDs)
marker_sets <- list(
  QC = c("AT3G11260"),  # WOX5
  Endodermis = c("AT3G54220", "AT4G37650", "AT5G57620", "AT5G62160"),  # SCR, SHR, MYB36, CASP1
  Cortex = c("AT5G14740", "AT3G17720"),  # CO2, PEP
  Trichoblast = c("AT5G49270", "AT1G69530", "AT5G12330", "AT4G25590"),  # COBL9, EXP7, RHD6, ADF8
  Atrichoblast = c("AT1G79840"),  # GL2
  Phloem = c("AT1G79430", "AT1G71890", "AT3G12810"),  # APL, SUC2, S17
  Xylem = c("AT1G71930", "AT4G35350", "AT5G05340"),  # VND7, XCP1, IRX1
  Procambium = c("AT4G32880", "AT2G01830"),  # ATHB8, WOL/CRE1
  Columella_LRC = c("AT1G62300", "AT5G67400")  # columella/LRC representatives
)
all_markers <- unique(unlist(marker_sets))
all_markers <- all_markers[all_markers %in% rownames(counts)]

seu <- CreateSeuratObject(counts = counts, project = "GSE123818_WT", min.cells = 3, min.features = 200)

# Replicate ID from 10x barcode suffix
seu$replicate <- paste0("rep", sub(".*-", "", colnames(seu)))
seu$SubjectName <- seu$replicate

# Organellar QC
org_genes <- rownames(seu)[grepl("^ATMG|^ATCG", rownames(seu))]
seu[["percent.org"]] <- PercentageFeatureSet(seu, features = org_genes)

# Filtering thresholds (stated in report)
# nFeature 500–8000, nCount >= 500, organellar <= 5%
seu <- subset(seu, subset = nFeature_RNA >= 500 & nFeature_RNA <= 8000 &
                nCount_RNA >= 500 & percent.org <= 5)

# Exclude organellar genes from HVG analysis
genes_use <- rownames(seu)[!grepl("^ATMG|^ATCG", rownames(seu))]
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000,
                            features = genes_use, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)

# Markers per cluster
Idents(seu) <- "seurat_clusters"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(ROOT, "tables", "sc_cluster_markers.csv"), row.names = FALSE)

# Score cell types and assign clusters
score_cols <- character()
for (nm in names(marker_sets)) {
  feats <- intersect(marker_sets[[nm]], rownames(seu))
  if (length(feats) >= 1) {
    col_nm <- paste0("score_", nm)
    seu <- AddModuleScore(seu, features = list(feats), name = col_nm, search = TRUE, verbose = FALSE)
    score_cols <- c(score_cols, paste0(col_nm, "1"))
  }
}

cluster_scores <- aggregate(
  seu@meta.data[, score_cols, drop = FALSE],
  by = list(cluster = seu$seurat_clusters),
  FUN = mean
)
rownames(cluster_scores) <- cluster_scores$cluster

best_type <- apply(cluster_scores[, score_cols, drop = FALSE], 1, function(x) {
  nm <- score_cols[which.max(x)]
  sub("^score_", "", sub("1$", "", nm))
})
cluster_map <- setNames(as.character(best_type), rownames(cluster_scores))
seu$cell_type <- unname(cluster_map[as.character(seu$seurat_clusters)])

write.csv(cluster_scores, file.path(ROOT, "tables", "sc_cluster_celltype_scores.csv"), row.names = FALSE)

# Plots
p1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(seu, group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle("Annotated cell types")
ggsave(file.path(ROOT, "figures", "sc_umap_clusters.png"), p1, width = 8, height = 6, dpi = 150)
ggsave(file.path(ROOT, "figures", "sc_umap_clusters.pdf"), p1, width = 8, height = 6)
ggsave(file.path(ROOT, "figures", "sc_umap_celltype.png"), p2, width = 9, height = 6, dpi = 150)
ggsave(file.path(ROOT, "figures", "sc_umap_celltype.pdf"), p2, width = 9, height = 6)

marker_present <- all_markers[all_markers %in% rownames(seu)]
if (length(marker_present) > 0) {
  p_dot <- DotPlot(seu, features = marker_present, group.by = "cell_type") +
    RotatedAxis() + ggtitle("Root marker genes")
  ggsave(file.path(ROOT, "figures", "sc_dotplot_markers.png"), p_dot, width = 12, height = 6, dpi = 150)
  ggsave(file.path(ROOT, "figures", "sc_dotplot_markers.pdf"), p_dot, width = 12, height = 6)

  p_vln <- VlnPlot(seu, features = head(marker_present, min(6, length(marker_present))),
                   group.by = "cell_type", pt.size = 0)
  ggsave(file.path(ROOT, "figures", "sc_violin_markers.png"), p_vln, width = 12, height = 8, dpi = 150)
  ggsave(file.path(ROOT, "figures", "sc_violin_markers.pdf"), p_vln, width = 12, height = 8)
}

saveRDS(seu, file.path(ROOT, "objects", "seurat_wt_annotated.rds"))

# Export for h5ad bridge (Matrix Market + metadata)
export_dir <- file.path(ROOT, "objects", "h5ad_export")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
counts_mat <- GetAssayData(seu, layer = "counts")
Matrix::writeMM(counts_mat, file.path(export_dir, "matrix.mtx"))
writeLines(colnames(counts_mat), file.path(export_dir, "barcodes.tsv"))
writeLines(rownames(counts_mat), file.path(export_dir, "features.tsv"))
write.csv(seu@meta.data, file.path(export_dir, "meta.csv"))
write.csv(Embeddings(seu, "umap"), file.path(export_dir, "umap.csv"))
write.csv(Embeddings(seu, "pca"), file.path(export_dir, "pca.csv"))
write.csv(data.frame(gene = VariableFeatures(seu)), file.path(export_dir, "hvg.csv"), row.names = FALSE)
data_mat <- GetAssayData(seu, layer = "data")
saveRDS(list(data = data_mat), file.path(ROOT, "objects", "seurat_data_layer.rds"))

message("Phase 1 complete. Cells: ", ncol(seu), " Clusters: ", length(unique(seu$seurat_clusters)))
