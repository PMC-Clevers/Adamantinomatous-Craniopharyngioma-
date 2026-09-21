################################################################################
## pseudotime_monocle3.R
## Monocle3 pseudotime on the dental-like ACP epithelium.
## Root = Transitioning; branch = Transitioning -> Ameloblast lineage.
## Stage-associated genes plotted over pseudotime (Figure 2G-I).
################################################################################

library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

dental_acp <- readRDS("dental_acp.rds")   # dental-like epithelial subset

## =============================================================================
## 1. Seurat -> cell_data_set, carrying the existing UMAP and clusters
## =============================================================================

cds <- as.cell_data_set(dental_acp)
cds <- estimate_size_factors(cds)                       # size-factor normalization
rowData(cds)$gene_short_name <- rownames(cds)           # needed for gene plots

## transfer Seurat UMAP + clustering into the cds
reducedDims(cds)$UMAP <- Embeddings(dental_acp, "umap.rpca")
cds@clusters$UMAP$clusters <- setNames(as.character(dental_acp$annotation),
                                       colnames(dental_acp))
cds@clusters$UMAP$partitions <- setNames(rep(1, ncol(cds)), colnames(cds))

## =============================================================================
## 2. Learn the principal graph and root in the Transitioning population
## =============================================================================

cds <- learn_graph(cds, use_partition = FALSE)

## pick the root principal node where Transitioning cells are most enriched
get_root_node <- function(cds, root_group, col = "annotation") {
  cell_ids <- which(colData(cds)[[col]] == root_group)
  vg <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  vertex_of_cell <- as.character(
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[, 1])
  names(vertex_of_cell) <- colnames(cds)
  tab <- table(vertex_of_cell[cell_ids])
  vg[as.numeric(names(which.max(tab)))]
}
root_node <- get_root_node(cds, "Transitioning")
cds <- order_cells(cds, root_pr_nodes = root_node)

## pseudotime + annotation UMAPs (Figure 2G, and coloured by cell type)
plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE,
           label_leaves = FALSE, label_roots = FALSE) + ggtitle("Pseudotime")
ggsave("pseudotime_UMAP.pdf", width = 6, height = 5)

plot_cells(cds, color_cells_by = "annotation", label_cell_groups = FALSE) +
  ggtitle("Annotation")
ggsave("pseudotime_annotation_UMAP.pdf", width = 7, height = 5)

## =============================================================================
## 3. Branch selection: Transitioning -> Ameloblast lineage
##    Interactive (choose the branch on the graph):
##      cds_branch <- choose_graph_segments(cds)
##    Reproducible alternative: restrict to the populations on this branch.
## =============================================================================

branch_labels <- c("Transitioning","Dental Epithelium","OEE/SR","IEE",
                   "Pre-ameloblast","Ameloblast")   # CONFIRM against your annotation
cds_branch <- cds[, colData(cds)$annotation %in% branch_labels]

## =============================================================================
## 4. Genes varying along the trajectory (graph_test)
## =============================================================================

graph_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
graph_res <- graph_res[order(-graph_res$morans_I), ]
write.csv(graph_res, "pseudotime_graphtest_genes.csv")

## =============================================================================
## 5. Stage-associated genes over pseudotime (Figure 2H-I)
##    Order kept as: early -> intermediate -> late (do not reorder).
## =============================================================================

stage_genes <- c("KRT15","PAPPA","SPINK5",     # early
                 "SOX6","SCUBE3","DLX5",        # intermediate
                 "SP6","KIF5C","VWDE","KRT17")  # late

plot_genes_in_pseudotime(
  cds_branch[stage_genes, ],
  color_cells_by = "annotation",
  min_expr = 0.5,
  ncol = 2
) + scale_color_manual(values = brewer.pal(max(3, length(branch_labels)), "Set2"))
ggsave("dental_Transitioning_to_Amelo_genes_over_pseudotime.pdf", width = 8, height = 12)

saveRDS(cds, "dental_acp_monocle3_cds.rds")
