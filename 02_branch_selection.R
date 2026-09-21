## =============================================================================
## 02_branch_selection.R
##
## Select a specific lineage branch (e.g. progenitor population -> terminal
## cell type) from a full cell_data_set trajectory, and root pseudotime
## within that branch specifically.
##
## choose_graph_segments() is INTERACTIVE - it opens a plot window and
## requires manually clicking along the branch of interest, then confirming
## the selection. This cannot be run non-interactively / in a batch script.
## =============================================================================

library(monocle3)
library(SingleCellExperiment)
library(ggplot2)

source("scripts/01_pseudotime_monocle3.R")   # for get_earliest_principal_node()

## -----------------------------------------------------------------------
## Generic branch-selection pipeline
##
## full_cds:      the full trajectory cell_data_set (e.g. dental_cds)
## root_population: population name to root the SELECTED BRANCH at (must
##                   itself be included in the manual selection)
## -----------------------------------------------------------------------
select_and_root_branch <- function(full_cds, root_population) {

  ## 1. Inspect the full graph first - identify where populations of
  ##    interest sit and whether a continuous path exists between them
  print(
    plot_cells(full_cds, color_cells_by = "annotation", label_cell_groups = TRUE,
               label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 3)
  )

  ## 2. Interactive selection - click along the path of interest, then "Done"
  cds_branch <- choose_graph_segments(full_cds)

  ## 3. Confirm the selection captured the expected populations
  message("Cells in selected branch: ", ncol(cds_branch))
  print(table(colData(cds_branch)$annotation))

  ## 4. Carry over UMAP coordinates (choose_graph_segments can drop this)
  reducedDim(cds_branch, "UMAP") <- reducedDim(full_cds, "UMAP")[colnames(cds_branch), ]

  ## 5. Re-cluster and re-learn the graph on just this subset
  cds_branch <- cluster_cells(cds_branch, reduction_method = "UMAP")
  cds_branch <- learn_graph(cds_branch, use_partition = FALSE)

  ## 6. Root at the specified population within this branch
  cds_branch <- order_cells(
    cds_branch,
    root_pr_nodes = get_earliest_principal_node(cds_branch, root_population)
  )

  ## 7. Sanity check pseudotime - Inf values indicate cells disconnected
  ##    from the root on the learned graph
  pt_summary <- summary(pseudotime(cds_branch))
  print(pt_summary)
  n_inf <- sum(!is.finite(pseudotime(cds_branch)))
  if (n_inf > 0) {
    warning(n_inf, " cells have non-finite pseudotime and should be filtered ",
            "before downstream analysis (e.g. GeneSwitches).")
  }

  cds_branch
}

## -----------------------------------------------------------------------
## Diagnostic: check for a library-depth confound before trusting a
## branch's downstream results. Large differences in nCount_RNA/nFeature_RNA
## between populations in a branch can produce switch-gene signatures
## dominated by housekeeping/ribosomal genes rather than genuine biology,
## especially in small/imbalanced branches.
## -----------------------------------------------------------------------
check_branch_depth_confound <- function(cds_branch) {
  library(dplyr)
  df <- as.data.frame(colData(cds_branch)) %>%
    group_by(annotation) %>%
    summarise(
      n = n(),
      median_nCount = median(nCount_RNA),
      median_nFeature = median(nFeature_RNA),
      .groups = "drop"
    )
  print(df)

  rho <- suppressWarnings(
    cor(colData(cds_branch)$nCount_RNA, pseudotime(cds_branch), method = "spearman")
  )
  message("Spearman correlation of nCount_RNA with pseudotime: ", round(rho, 3))
  if (!is.na(rho) && abs(rho) > 0.3) {
    warning("Moderate-to-strong correlation between library depth and pseudotime ",
            "detected. Switch-gene results from this branch may be confounded by ",
            "sequencing depth rather than reflecting genuine differentiation ",
            "biology - interpret with caution, particularly if population sizes ",
            "are also highly imbalanced.")
  }
  invisible(df)
}

## =============================================================================
## Example usage: Dental, Dental Epithelium -> Ameloblast branch
## =============================================================================

dental_cds <- readRDS("dental_cds_monocle3.rds")

cds_DE_to_Amelo <- select_and_root_branch(dental_cds, root_population = "Transitioning")
check_branch_depth_confound(cds_TR_to_Amelo)

plot_cells(cds_TR_to_Amelo, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3)
ggsave("output/dental_TR_to_Amelo_pseudotime.pdf", width = 7, height = 6)

saveRDS(cds_DE_to_Amelo, "dental_TR_to_Amelo_cds.rds")

## Repeat for other branches of interest, e.g.:
## cds_Transitioning_to_Amelo <- select_and_root_branch(dental_cds, "Dental epithelium")
## cds_CervicalLoop_to_Amelo  <- select_and_root_branch(dental_cds, "Cervical Loop")
##
## Give each branch object a distinct name and save it immediately -
## do not reuse a single variable name (e.g. cds_subset) across multiple
## branch selections, as re-rooting or reselecting will silently overwrite
## the previous branch's result.
