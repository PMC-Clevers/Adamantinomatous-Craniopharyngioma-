## =============================================================================
## 01_pseudotime_monocle3.R
##
## Build Monocle3 cell_data_set objects from annotated Seurat objects and
## order cells in pseudotime from a specified root population.
##
## Input:  Seurat objects with a UMAP reduction and an `annotation` metadata
##         column (e.g. oral_acp, dental_acp)
## Output: cell_data_set objects with pseudotime assigned, saved to disk
## =============================================================================

library(monocle3)
library(SeuratWrappers)
library(SingleCellExperiment)
library(ggplot2)

## -----------------------------------------------------------------------
## Helper: identify the principal graph node closest to a given root
## population, used to seed order_cells()
## -----------------------------------------------------------------------
get_earliest_principal_node <- function(cds, root_cluster) {
  cell_ids <- which(colData(cds)$annotation == root_cluster)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  root_pr_nodes
}

## -----------------------------------------------------------------------
## Generic pipeline: Seurat object -> rooted, ordered cell_data_set
##
## reduction_slot_pattern: regex matching the UMAP reduction name in the
##   Seurat object (varies between objects, e.g. "UMAP.DIM20" vs "umap1") -
##   inspect reducedDimNames(cds) after the initial conversion if unsure.
## -----------------------------------------------------------------------
build_and_order_cds <- function(seurat_obj, root_population,
                                 reduction_slot_pattern = "^UMAP") {

  cds <- as.cell_data_set(seurat_obj)
  cds@clusters$UMAP$clusters <- seurat_obj$annotation

  ## normalize whatever the UMAP reduction is named to "UMAP"
  rd_names <- names(cds@int_colData@listData$reducedDims)
  match_idx <- grepl(reduction_slot_pattern, rd_names)
  if (!any(match_idx)) {
    stop("No reducedDim name matched pattern '", reduction_slot_pattern,
         "'. Available: ", paste(rd_names, collapse = ", "))
  }
  names(cds@int_colData@listData$reducedDims)[match_idx] <- "UMAP"

  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = FALSE)

  cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, root_population))

  cds
}

## =============================================================================
## Example usage - Oral
## =============================================================================

oral_acp <- readRDS("oral_acp_with_pseudotime_early_oral_root.rds")

oral_cds <- build_and_order_cds(
  oral_acp,
  root_population = "Early-Oral epithelium",
  reduction_slot_pattern = "^UMAP"   # e.g. matches "UMAP.DIM20"
)

plot_cells(oral_cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3)
ggsave("output/pseudotime_oral_early_oral_root.pdf", width = 7, height = 6)

oral_acp$pseudotime <- pseudotime(oral_cds)
saveRDS(oral_acp, "oral_acp_with_pseudotime_early_oral_root.rds")
saveRDS(oral_cds, "oral_cds_monocle3.rds")

## =============================================================================
## Example usage - Dental
## =============================================================================

dental_acp <- readRDS("dental_acp_with_pseudotime_transitioning_root.rds")

## check available root populations first
table(dental_acp$annotation)

dental_cds <- build_and_order_cds(
  dental_acp,
  root_population = "Transitioning",
  reduction_slot_pattern = "^UMAP"
)

plot_cells(dental_cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3)
ggsave("output/pseudotime_dental_transitioning_root.pdf", width = 7, height = 6)

dental_acp$pseudotime <- pseudotime(dental_cds)
saveRDS(dental_acp, "dental_acp_with_pseudotime_transitioning_root.rds")
saveRDS(dental_cds, "dental_cds_monocle3.rds")

## -----------------------------------------------------------------------
## NOTE on re-rooting an existing branch object
## -----------------------------------------------------------------------
## order_cells() recomputes pseudotime as geodesic distance along the
## EXISTING principal graph from a new root - it does not rebuild the graph
## itself. If a cds has been re-rooted multiple times across a session,
## rebuild it fresh (cluster_cells() + learn_graph()) from a saved,
## known-good state before re-rooting again, rather than repeatedly
## re-rooting the same live object. Confirm the result with:
##   summary(pseudotime(cds))
##   plot_cells(cds, color_cells_by = "pseudotime", ...)
## before proceeding to downstream analysis.

## -----------------------------------------------------------------------
## NOTE on annotation changes during revision
## -----------------------------------------------------------------------
## If a population's annotation is refined or relabeled during revision
## (e.g. based on additional marker validation), document:
##   (1) the specific markers/evidence supporting the change
##   (2) that the change is applied consistently across all downstream
##       figures/tables, not just the one a reviewer commented on
## Do not choose a root population based on which choice produces a more
## favorable-looking trajectory result.
