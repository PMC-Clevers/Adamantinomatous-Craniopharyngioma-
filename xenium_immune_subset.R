################################################################################
## Xenium immune subset - subclustering and annotation (annotation_immune)
## Produces xenium_with_immune_subtypes.rds, used for the immune spatial analysis.
## Runs on the integrated Xenium object from xenium_analysis_final.R.
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(scCustomize)
library(dplyr)
library(future)
options(future.globals.maxSize = 60 * 1024^3)

xenium <- readRDS("20251103_xenium-nobt259.rds")

################################################################################
## 1. Subset the immune / myeloid compartment
##    CONFIRM: which parent clusters were immune. In annotation_v2 the immune
##    cells sit in the "Myeloid" clusters; T cells are pulled out by reclustering.
##    If you also had a lymphoid parent cluster, add it here.
################################################################################

immune <- subset(xenium, subset = annotation_v2 %in% c("Myeloid"))
# (alternative if you subset on PTPRC+ instead:
#  immune <- subset(xenium, subset = PTPRC > 0) )

################################################################################
## 2. Re-normalize, integrate and cluster the immune subset
################################################################################

DefaultAssay(immune) <- "Xenium"
immune <- SCTransform(immune, assay = "Xenium", conserve.memory = TRUE, verbose = FALSE)
immune <- RunPCA(immune, npcs = 30, verbose = FALSE)
immune <- RunHarmony(immune, group.by.vars = "sample", assay.use = "SCT",
                     reduction = "pca", dims.use = 1:30, reduction.save = "harmony")
immune <- RunUMAP(immune, reduction = "harmony", dims = 1:20)
immune <- FindNeighbors(immune, reduction = "harmony", dims = 1:20)
immune <- FindClusters(immune, resolution = 0.5)   # CONFIRM resolution used

################################################################################
## 3. Markers used for annotation
################################################################################

immune_markers <- c(
  # pan-immune / macrophage
  "PTPRC","CD68","CSF1R","LYZ","CD163","MRC1","APOE","SELENOP","TMEM176A","TMEM176B",
  # GPNMB+/APOE+ lipid-associated TAM
  "GPNMB","PLIN2",
  # M1 / inflammatory
  "CD80","CD86","IL1B","CXCL8","STAT1","IRF1","GBP2",
  # CCL2+ / monocytic
  "CCL2","S100A8","S100A9","VCAN","CCR2",
  # T cells
  "CD3D","CD3E","CD8A","CD4","CCR7","FOXP3",
  # B / plasma
  "CD79A","MS4A1","CD38","SDC1")
immune_markers <- immune_markers[immune_markers %in% rownames(immune)]

## inspect cluster markers to assign labels
DefaultAssay(immune) <- "SCT"
for (i in seq_along(immune[["SCT"]]@SCTModel.list))
  immune[["SCT"]]@SCTModel.list[[i]]@umi.assay <- "Xenium"
future::plan("sequential")
immune <- PrepSCTFindMarkers(immune, assay = "SCT", verbose = TRUE)
imm_markers <- FindAllMarkers(immune, assay = "SCT", only.pos = TRUE,
                              min.pct = 0.1, logfc.threshold = 0.5)
write.csv(imm_markers, "xenium_immune_markers.csv")

################################################################################
## 4. Annotate immune clusters (annotation_immune)
##    CONFIRM the cluster -> label map against your dotplot; the labels used
##    in the paper are the six below (+ Low-quality, dropped for spatial).
################################################################################

immune_ids <- c(
  # example map, EDIT numbers to match your clusters:
  "0" = "T-cells (CD8)",
  "1" = "Macrophages",
  "2" = "TAM (GPNMB+/APOE+)",
  "3" = "T-cells (CD4/Treg)",
  "4" = "Macrophages (CCL2+)",
  "5" = "Macrophages (CD86+)",
  "6" = "Low-quality"
)
immune$annotation_immune <- plyr::mapvalues(as.character(Idents(immune)),
                                            from = names(immune_ids),
                                            to = immune_ids, warn_missing = FALSE)
Idents(immune) <- "annotation_immune"
saveRDS(immune, "xenium_myeloid_annotated.rds")

################################################################################
## 5. UMAP + marker dot plot (Supplementary Figure S4/S10A-B)
################################################################################

DimPlot_scCustom(immune, group.by = "annotation_immune") +
  ggtitle("Xenium immune subset")
ggsave("xenium_immune_UMAP.pdf", width = 8, height = 6)

DotPlot_scCustom(immune, features = immune_markers, group.by = "annotation_immune") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Xenium immune subset annotation")
ggsave("xenium_immune_dotplot.pdf", width = 12, height = 5)

################################################################################
## 6. Merge immune labels back into the full object (annotation_immune)
##    Non-immune cells keep their annotation_v2 label; immune cells get the
##    fine subtype. This is the object used for the spatial (Squidpy) analysis.
################################################################################

lab <- setNames(as.character(immune$annotation_immune), colnames(immune))
xenium$annotation_immune <- as.character(xenium$annotation_v2)
xenium$annotation_immune[names(lab)] <- lab

saveRDS(xenium, "xenium_with_immune_subtypes.rds")
