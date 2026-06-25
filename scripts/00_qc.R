#!/usr/bin/env Rscript
# Phase 0 — Setup & input QC
set.seed(42)
ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
dir.create(file.path(ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "objects"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "logs"), showWarnings = FALSE, recursive = TRUE)

read_gz_csv <- function(path, sep = ",", gene_col = 1) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  df <- read.table(con, sep = sep, header = TRUE, check.names = FALSE, quote = "", comment.char = "")
  rownames(df) <- df[[gene_col]]
  df <- df[, -gene_col, drop = FALSE]
  as.data.frame(df)
}

sc_path <- file.path(ROOT, "data", "GSE123818_Root_single_cell_wt_datamatrix.csv.gz")
bulk_path <- file.path(ROOT, "data", "GSE123818_Root_bulk_tissue_datamatrix.txt.gz")

message("Loading single-cell matrix...")
sc <- read_gz_csv(sc_path)

message("Loading bulk matrix...")
bulk <- read_gz_csv(bulk_path, sep = "\t", gene_col = 1)

count_cols <- grep("\\.readcount$", colnames(bulk), value = TRUE)
bulk_counts <- bulk[, count_cols, drop = FALSE]
colnames(bulk_counts) <- gsub("\\.readcount$", "", colnames(bulk_counts))

shared_genes <- intersect(rownames(sc), rownames(bulk_counts))

sc_sparsity <- mean(as.matrix(sc) == 0) * 100
bulk_sparsity <- mean(as.matrix(bulk_counts) == 0) * 100

qc_summary <- data.frame(
  dataset = c("single_cell_wt", "bulk_tissue"),
  n_genes = c(nrow(sc), nrow(bulk_counts)),
  n_features = c(ncol(sc), ncol(bulk_counts)),
  gene_id_example = c(rownames(sc)[1], rownames(bulk_counts)[1]),
  pct_zero = c(sc_sparsity, bulk_sparsity),
  n_shared_genes = c(length(shared_genes), length(shared_genes)),
  stringsAsFactors = FALSE
)

write.csv(qc_summary, file.path(ROOT, "tables", "phase0_qc_summary.csv"), row.names = FALSE)
write.csv(data.frame(gene = shared_genes), file.path(ROOT, "tables", "shared_genes.csv"), row.names = FALSE)

saveRDS(list(sc = sc, bulk = bulk_counts, bulk_annot = bulk, shared_genes = shared_genes),
        file.path(ROOT, "objects", "phase0_input.rds"))

sink(file.path(ROOT, "logs", "sessionInfo_phase0.txt"))
print(sessionInfo())
sink()

message("Phase 0 complete.")
print(qc_summary)
