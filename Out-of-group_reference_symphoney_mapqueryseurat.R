## ============================================================
## FINAL SCRIPT: Out-of-group reference mapping
## Seurat MapQuery + Symphony, 4-tissue combined reference
## (Dental, Gingival, Epidermis, Esophagus - Tongue excluded)
## ============================================================

library(Seurat)
library(symphony)
library(harmony)      # must be v0.1.1 for Symphony compatibility
library(irlba)
library(uwot)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

setwd("E:/Laurens")

tissue_colors <- c(Dental = "#F8766D", Epidermis = "#7CAE00", Esophagus = "#00BFC4", Gingival = "#C77CFF")

## ============================================================
## PART 0: Load your 5 base objects (adjust paths as needed)
## Each must already be: SCTransform'd, RunPCA'd, RunUMAP'd (return.model=TRUE)
## and have a working RNA "data" layer (NormalizeData if empty)
## ============================================================

dental_acp                    <- readRDS("Dental_subset_tissue_scRNAseq_SCT_PCA_UMAPrpca2_clusters.rds")
oral_acp                      <- readRDS("Oral_subset_tissue_scRNAseq_SCT_PCA_UMAPrpca2_clusters.rds")
dental_epithelial_lineage_so  <- readRDS("dental_epithelial_lineage_so_PCA_FIXED.rds")
gingival_epi                  <- readRDS("gingival_epi_prepped_FINAL.rds")
skin_epi                      <- readRDS("skin_epi_map.rds")
esophagus_epi                 <- readRDS("esophagus_epi.rds")

## ensure RNA "data" layer is populated on ACP query objects (fixes "Layer data is empty")
if (nrow(GetAssayData(dental_acp, assay = "RNA", layer = "data")) == 0) {
  dental_acp <- NormalizeData(dental_acp, assay = "RNA")
}
if (nrow(GetAssayData(oral_acp, assay = "RNA", layer = "data")) == 0) {
  oral_acp <- NormalizeData(oral_acp, assay = "RNA")
}

## tag tissue identity on each reference
dental_epithelial_lineage_so$tissue <- "Dental"
gingival_epi$tissue                 <- "Gingival"
skin_epi$tissue                     <- "Epidermis"
esophagus_epi$tissue                <- "Esophagus"

## ensure consistent assay naming (fix RNA_sym -> RNA if present)
fix_assay_name <- function(obj) {
  if ("RNA_sym" %in% Assays(obj)) {
    obj[["RNA"]] <- obj[["RNA_sym"]]
    obj[["RNA_sym"]] <- NULL
    DefaultAssay(obj) <- "RNA"
  }
  obj
}
dental_epithelial_lineage_so <- fix_assay_name(dental_epithelial_lineage_so)
gingival_epi                 <- fix_assay_name(gingival_epi)
skin_epi                     <- fix_assay_name(skin_epi)
esophagus_epi                <- fix_assay_name(esophagus_epi)


## ============================================================
## PART 1: SEURAT MAPQUERY
## ============================================================

## ---- 1a. Build combined 4-tissue reference (Seurat) ----
combined_ref_no_tongue <- merge(dental_epithelial_lineage_so,
                                y = list(gingival_epi, skin_epi, esophagus_epi))

DefaultAssay(combined_ref_no_tongue) <- "RNA"
if ("SCT" %in% Assays(combined_ref_no_tongue)) combined_ref_no_tongue[["SCT"]] <- NULL
combined_ref_no_tongue <- JoinLayers(combined_ref_no_tongue)

ncol(combined_ref_no_tongue)   # should be 149,258
table(combined_ref_no_tongue$tissue)

combined_ref_no_tongue <- SCTransform(combined_ref_no_tongue, conserve.memory = TRUE, verbose = TRUE)
combined_ref_no_tongue <- RunPCA(combined_ref_no_tongue, npcs = 30, verbose = TRUE)
combined_ref_no_tongue <- RunUMAP(combined_ref_no_tongue, reduction = "pca", dims = 1:30, return.model = TRUE)

saveRDS(combined_ref_no_tongue, "combined_ref_no_tongue_SCT.rds")

## ---- 1b. map_and_score function (unintegrated PCA-based; avoids SCT-model version bug) ----
map_and_score <- function(reference, query, ref_name, query_name, label_col) {
  shared_feats <- intersect(VariableFeatures(reference), VariableFeatures(query))
  
  anchors <- FindTransferAnchors(
    reference            = reference,
    query                = query,
    normalization.method = "SCT",
    reference.assay      = "SCT",
    query.assay          = "SCT",
    reference.reduction  = "pca",
    dims                 = 1:30,
    features             = shared_feats,
    verbose              = TRUE
  )
  
  mapped <- MapQuery(
    anchorset           = anchors,
    query               = query,
    reference           = reference,
    reference.reduction = "pca",
    reduction.model      = "umap",
    refdata              = list(predicted.celltype = label_col)
  )
  
  mapped <- AddMetaData(mapped, MappingScore(anchors, ndim = 30), col.name = "mapping.score")
  mapped$reference_used <- ref_name
  mapped$query_used     <- query_name
  mapped$n_shared_feats <- length(shared_feats)
  mapped
}

## ---- 1c. Run mapping ----
mapped_dental_notongue_seurat <- map_and_score(combined_ref_no_tongue, dental_acp, "Combined4Tissue", "Dental_ACP", label_col = "tissue")
saveRDS(mapped_dental_notongue_seurat, "seurat_mapquery_dental_4tissue_notongue.rds")

mapped_oral_notongue_seurat <- map_and_score(combined_ref_no_tongue, oral_acp, "Combined4Tissue", "Oral_ACP", label_col = "tissue")
saveRDS(mapped_oral_notongue_seurat, "seurat_mapquery_oral_4tissue_notongue.rds")

table(mapped_dental_notongue_seurat$predicted.predicted.celltype)
table(mapped_oral_notongue_seurat$predicted.predicted.celltype)


## ============================================================
## PART 2: SYMPHONY
## ============================================================

## ---- 2a. Build combined 4-tissue reference (Symphony) ----
combined_ref_symphony <- merge(dental_epithelial_lineage_so,
                               y = list(gingival_epi, skin_epi, esophagus_epi))
DefaultAssay(combined_ref_symphony) <- "RNA"
if ("SCT" %in% Assays(combined_ref_symphony)) combined_ref_symphony[["SCT"]] <- NULL
combined_ref_symphony <- JoinLayers(combined_ref_symphony)
combined_ref_symphony <- NormalizeData(combined_ref_symphony, assay = "RNA")

expr_mat  <- GetAssayData(combined_ref_symphony, assay = "RNA", layer = "data")
meta_data <- combined_ref_symphony@meta.data
stopifnot(identical(colnames(expr_mat), rownames(meta_data)))

vargenes_res <- symphony::vargenes_vst(expr_mat, groups = rep("all", ncol(expr_mat)), topn = 2000)
gene_means <- Matrix::rowMeans(expr_mat[vargenes_res, ])
gene_sds   <- symphony::rowSDs(expr_mat[vargenes_res, ], gene_means)
vargenes_means_sds <- tibble(symbol = vargenes_res, mean = gene_means, stddev = gene_sds)

exp_ref_scaled <- symphony::scaleDataWithStats(expr_mat[vargenes_res, ], gene_means, gene_sds, 1)

set.seed(0)
s <- irlba::irlba(exp_ref_scaled, nv = 30)
Z_pca_ref <- diag(s$d) %*% t(s$v)
loadings  <- s$u

ref_harmObj <- harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref), meta_data = meta_data, theta = c(2),
  vars_use = c('tissue'), nclust = 100, max.iter.harmony = 20,
  return_object = TRUE, do_pca = FALSE
)

## build Symphony reference manually (bypasses getZorig/getZcorr version bug)
set.seed(111)
res <- list(meta_data = meta_data, vargenes = vargenes_means_sds, loadings = loadings)
res$R <- ref_harmObj$R
res$Z_orig <- ref_harmObj$Z_orig
res$Z_corr <- ref_harmObj$Z_corr
res$betas  <- ref_harmObj$W
res$centroids <- t(symphony:::cosine_normalize_cpp(ref_harmObj$R %*% t(ref_harmObj$Z_corr), 1))
res$cache <- symphony:::compute_ref_cache(res$R, res$Z_corr)

colnames(res$Z_orig) <- row.names(meta_data)
rownames(res$Z_orig) <- paste0("PC_", seq_len(nrow(res$Z_corr)))
colnames(res$Z_corr) <- row.names(meta_data)
rownames(res$Z_corr) <- paste0("harmony_", seq_len(nrow(res$Z_corr)))

cluster_sizes <- res$cache[[1]] %>% as.matrix()
centroid_sums <- t(res$Z_corr %*% t(res$R)) %>% as.data.frame()
centroids_pc  <- sweep(centroid_sums, 1, cluster_sizes, "/")
colnames(centroids_pc) <- paste0("harmony_", seq_len(nrow(res$Z_corr)))
rownames(centroids_pc) <- paste0("centroid_", seq_len(nrow(res$R)))
res$centroids_pc <- centroids_pc

umap <- uwot::umap(t(res$Z_corr), n_neighbors = 30, learning_rate = 0.5,
                   init = "laplacian", metric = "cosine", fast_sgd = FALSE,
                   n_sgd_threads = 1, min_dist = 0.1, n_threads = 4, ret_model = TRUE)
res$umap$embedding <- umap$embedding
colnames(res$umap$embedding) <- c("UMAP1", "UMAP2")

if (file.exists('E:/Laurens/symphony_uwot_model_notongue')) file.remove('E:/Laurens/symphony_uwot_model_notongue')
uwot::save_uwot(umap, file = 'E:/Laurens/symphony_uwot_model_notongue', unload = FALSE, verbose = FALSE)
res$save_uwot_path <- 'E:/Laurens/symphony_uwot_model_notongue'

symphony_ref_notongue <- res
saveRDS(symphony_ref_notongue, "symphony_combined_ref_notongue.rds")

## ---- 2b. Map ACP subsets onto Symphony reference ----
map_symphony <- function(query_obj, ref_obj, k = 5) {
  q_expr <- GetAssayData(query_obj, assay = "RNA", layer = "data")
  q_meta <- query_obj@meta.data
  mapped <- symphony::mapQuery(q_expr, q_meta, ref_obj, vars = NULL, do_normalize = FALSE)
  mapped <- symphony::knnPredict(mapped, ref_obj, ref_obj$meta_data$tissue, k = k)
  mapped
}

mapped_dental_4tissue <- map_symphony(dental_acp, symphony_ref_notongue)
saveRDS(mapped_dental_4tissue, "symphony_dental_acp_mapped_notongue.rds")

mapped_oral_4tissue <- map_symphony(oral_acp, symphony_ref_notongue)
saveRDS(mapped_oral_4tissue, "symphony_oral_acp_mapped_notongue.rds")

table(mapped_dental_4tissue$meta_data$cell_type_pred_knn)
table(mapped_oral_4tissue$meta_data$cell_type_pred_knn)


## ============================================================
## PART 3: SUMMARY TABLE + PLOTS
## ============================================================

build_table <- function(pred_vector, query_name, method_name) {
  df <- as.data.frame(table(pred_vector))
  colnames(df) <- c("predicted_tissue", "n")
  df$proportion <- round(df$n / sum(df$n) * 100, 1)
  df$query <- query_name
  df$method <- method_name
  df
}

master_table <- bind_rows(
  build_table(mapped_dental_notongue_seurat$predicted.predicted.celltype, "Dental_ACP", "Seurat MapQuery"),
  build_table(mapped_oral_notongue_seurat$predicted.predicted.celltype, "Oral_ACP", "Seurat MapQuery"),
  build_table(mapped_dental_4tissue$meta_data$cell_type_pred_knn, "Dental_ACP", "Symphony"),
  build_table(mapped_oral_4tissue$meta_data$cell_type_pred_knn, "Oral_ACP", "Symphony")
) %>%
  select(method, query, predicted_tissue, n, proportion) %>%
  arrange(method, query, desc(n))

print(master_table)
write.csv(master_table, "master_prediction_table_both_methods.csv", row.names = FALSE)

## Panel A - reference UMAPs
ref_df_seurat <- data.frame(Embeddings(combined_ref_no_tongue, "umap"), tissue = combined_ref_no_tongue$tissue)
colnames(ref_df_seurat)[1:2] <- c("UMAP_1", "UMAP_2")
p_ref_seurat <- ggplot(ref_df_seurat, aes(UMAP_1, UMAP_2, color = tissue)) +
  geom_point(size = 0.4, alpha = 0.7) + scale_color_manual(values = tissue_colors) +
  theme_classic() + ggtitle("Seurat: unintegrated combined reference")

ref_df_symphony <- data.frame(symphony_ref_notongue$umap$embedding, tissue = symphony_ref_notongue$meta_data$tissue)
colnames(ref_df_symphony)[1:2] <- c("UMAP_1", "UMAP_2")
p_ref_symphony <- ggplot(ref_df_symphony, aes(UMAP_1, UMAP_2, color = tissue)) +
  geom_point(size = 0.4, alpha = 0.7) + scale_color_manual(values = tissue_colors) +
  theme_classic() + ggtitle("Symphony: Harmony-integrated combined reference")

## Panel B/C - query overlays
qry_dental_seurat <- data.frame(Embeddings(mapped_dental_notongue_seurat, "ref.umap")); colnames(qry_dental_seurat)[1:2] <- c("UMAP_1","UMAP_2")
qry_oral_seurat   <- data.frame(Embeddings(mapped_oral_notongue_seurat, "ref.umap"));   colnames(qry_oral_seurat)[1:2]   <- c("UMAP_1","UMAP_2")
qry_dental_symphony <- data.frame(mapped_dental_4tissue$umap); colnames(qry_dental_symphony)[1:2] <- c("UMAP_1","UMAP_2")
qry_oral_symphony   <- data.frame(mapped_oral_4tissue$umap);   colnames(qry_oral_symphony)[1:2]   <- c("UMAP_1","UMAP_2")

p1 <- ggplot() + geom_point(data=ref_df_seurat, aes(UMAP_1,UMAP_2), color="gray85", size=0.3, alpha=0.4) +
  geom_point(data=qry_dental_seurat, aes(UMAP_1,UMAP_2), color="firebrick", size=0.5) + theme_classic() + ggtitle("Seurat: Dental_ACP")
p2 <- ggplot() + geom_point(data=ref_df_seurat, aes(UMAP_1,UMAP_2), color="gray85", size=0.3, alpha=0.4) +
  geom_point(data=qry_oral_seurat, aes(UMAP_1,UMAP_2), color="steelblue", size=0.5) + theme_classic() + ggtitle("Seurat: Oral_ACP")
p3 <- ggplot() + geom_point(data=ref_df_symphony, aes(UMAP_1,UMAP_2), color="gray85", size=0.3, alpha=0.4) +
  geom_point(data=qry_dental_symphony, aes(UMAP_1,UMAP_2), color="firebrick", size=0.5) + theme_classic() + ggtitle("Symphony: Dental_ACP")
p4 <- ggplot() + geom_point(data=ref_df_symphony, aes(UMAP_1,UMAP_2), color="gray85", size=0.3, alpha=0.4) +
  geom_point(data=qry_oral_symphony, aes(UMAP_1,UMAP_2), color="steelblue", size=0.5) + theme_classic() + ggtitle("Symphony: Oral_ACP")

(p_ref_seurat | p_ref_symphony) / (p1 | p2) / (p3 | p4)
ggsave("final_figure_all_panels.pdf", width = 12, height = 15)

## Panel D - stacked bar
ggplot(master_table, aes(x = query, y = proportion, fill = predicted_tissue)) +
  geom_col() +
  geom_text(aes(label = paste0(proportion, "%")), position = position_stack(vjust = 0.5), size = 5, color = "white") +
  facet_wrap(~method) +
  scale_fill_manual(values = tissue_colors) +
  theme_classic() +
  labs(title = "Predicted tissue composition: Seurat vs. Symphony (4-tissue reference)",
       x = NULL, y = "Proportion of cells (%)", fill = "Predicted tissue") +
  theme(plot.title = element_text(size = 18, face = "bold"), strip.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14))
ggsave("master_prediction_barplot_both_methods.pdf", width = 11, height = 7)