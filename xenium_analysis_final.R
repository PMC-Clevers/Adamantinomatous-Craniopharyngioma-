################################################################################
## Xenium analysis - final pipeline (paper version)
## Adamantinomatous craniopharyngioma (ACP)
## Only the steps that produced the objects and figures used in the manuscript.
## Runs kept: TMA + BT288 (BT259 excluded). QC: nCount_Xenium > 10.
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)          # for mapvalues (load before dplyr-heavy steps if needed)
library(RColorBrewer)
library(scCustomize)
library(future)
options(future.globals.maxSize = 60 * 1024^3)

setwd("E:/Laurens/Xenium")

################################################################################
## 1. Load runs and QC filter (nCount_Xenium > 10)
################################################################################

xenium.obj  <- LoadXenium("./output-XETG00101__0053865__Region_1__20251016__110509/",
                          fov = "fov", segmentations = "cell", flip.xy = TRUE)   # TMA
xenium.obj2 <- LoadXenium("./output-XETG00101__0053877__Region_1__20251016__110509/",
                          fov = "fov", segmentations = "cell", flip.xy = TRUE)   # BT288

## --- QC cell counts for Methods (before vs after) ---
qc_counts <- function(obj, name, min_counts = 10) {
  nb <- ncol(obj); na <- sum(obj$nCount_Xenium > min_counts)
  data.frame(run = name, before = nb, after = na,
             removed = nb - na, pct_removed = round(100 * (nb - na) / nb, 2))
}
qc <- rbind(qc_counts(xenium.obj, "TMA"), qc_counts(xenium.obj2, "BT288"))
qc <- rbind(qc, data.frame(run = "TOTAL", before = sum(qc$before), after = sum(qc$after),
                           removed = sum(qc$removed),
                           pct_removed = round(100 * sum(qc$removed) / sum(qc$before), 2)))
print(qc)

## --- apply filter ---
xenium.obj  <- subset(xenium.obj,  subset = nCount_Xenium > 10)
xenium.obj2 <- subset(xenium.obj2, subset = nCount_Xenium > 10)

################################################################################
## 2. Merge, normalize, integrate (Harmony by run) and cluster
################################################################################

xenium <- merge(xenium.obj, y = xenium.obj2)

## sample label from barcode prefix
xenium$sample <- sub("_.*$", "", colnames(xenium))
print(table(xenium$sample))

DefaultAssay(xenium) <- "Xenium"
xenium <- SCTransform(xenium, assay = "Xenium", conserve.memory = TRUE, verbose = FALSE)
xenium <- RunPCA(xenium, npcs = 50, verbose = FALSE)

xenium <- RunHarmony(xenium, group.by.vars = "sample", assay.use = "SCT",
                     reduction = "pca", dims.use = 1:50, theta = 4,
                     reduction.save = "harmony")

xenium <- RunUMAP(xenium, reduction = "harmony", dims = 1:25)
xenium <- FindNeighbors(xenium, reduction = "harmony", dims = 1:25)
xenium <- FindClusters(xenium, resolution = 0.6)   # -> SCT_snn_res.0.6, used for annotation_v2

################################################################################
## 3. Marker detection and annotation (annotation_v2)
################################################################################

## align SCT model UMI assay names before PrepSCTFindMarkers (v5 multi-model)
DefaultAssay(xenium) <- "SCT"
for (i in seq_along(xenium[["SCT"]]@SCTModel.list))
  xenium[["SCT"]]@SCTModel.list[[i]]@umi.assay <- "Xenium"
future::plan("sequential")
xenium <- PrepSCTFindMarkers(xenium, assay = "SCT", verbose = TRUE)

markers <- FindAllMarkers(xenium, assay = "SCT", only.pos = TRUE,
                          min.pct = 0.05, logfc.threshold = 0.5)
write.csv(markers, "xenium_markers_res06.csv")

## cluster -> label map (res 0.6, 29 clusters)
cluster.ids_v2 <- c(
  "Oral Epithelium 1", "Oral Epithelium 2", "Dental Epithelium / OEE", "Oral Epithelium 3",
  "Transitioning", "SI / SR", "Dental Epithelium", "Junctional basal", "Sulcus basal",
  "ATF3", "Ameloblast", "Glial", "PDL", "Pre-ameloblast", "Endothelial", "Myeloid",
  "Cervical Loop", "LYZ oral", "Myeloid", "IGFBP5", "Oral basal", "Junctional Keratinocytes",
  "Stress", "Sulcus Keratinocytes", "Cluster Cells", "CCL2 inflamed keratinocytes",
  "SRs", "Cycling", "BEST3 / SOX2")
names(cluster.ids_v2) <- as.character(0:28)

xenium$annotation_v2 <- plyr::mapvalues(as.character(Idents(xenium)),
                                        from = names(cluster.ids_v2),
                                        to = cluster.ids_v2, warn_missing = FALSE)
Idents(xenium) <- "annotation_v2"

saveRDS(xenium, "20251103_xenium-nobt259.rds")

## 4. Split ROIs into separate objects (per cropped image)
## coordinates were selected in xenium explorer and used for Cropping according to seurat instructions
################################################################################

outdir <- "Split_ROIs2"; dir.create(outdir, showWarnings = FALSE)
for (nm in names(xenium@images)) {
  ids <- tryCatch(xenium@images[[nm]]$centroids@cells, error = function(e) NULL)
  if (is.null(ids) || !length(ids)) next
  ids <- intersect(ids, colnames(xenium))
  so <- subset(xenium, cells = ids)
  so@images <- so@images[nm, drop = FALSE]
  saveRDS(so, file.path(outdir, paste0(nm, ".rds")))
}

## ROI objects used in the paper: all 15 (13 TMA cores + 2 BT288 pieces).
## TMA5 was not cropped (skipped), leaving TMA1-4 and TMA6-14 = 13 cores,
## plus ROI_BT288_1 and ROI_20251030_newregions_bt288_coordinates.
roi_files <- list.files("Split_ROIs2", pattern = "\\.rds$", full.names = TRUE)
roi_list  <- setNames(lapply(roi_files, readRDS),
                      sub("\\.rds$", "", basename(roi_files)))

length(roi_list)   # should be 15
names(roi_list)    # check ROI names
# reference any ROI as roi_list[["ROI_TMA11"]], roi_list[["ROI_BT288_1"]], etc.

################################################################################
## 5. Colour palette (stable label -> colour map)
################################################################################

parade_cols <- c(
  "#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000",
  "#7E6148","#B09C85","#FFDC91","#EFC000","#7f7fff","#af8dc3","#f1a340","#998ec3",
  "#fdb863","#b57259","#ff9da7","#5c7e19","#3c91e6","#f25f5c","#00a6a6","#ffb627",
  "#b8336a","cyan","purple","green","yellow")
labs_all   <- unique(cluster.ids_v2)
parade_map <- setNames(colorRampPalette(parade_cols)(length(labs_all)), labs_all)

################################################################################
## 6. Spatial cluster maps (Figure 6 / ImageDimPlots)
################################################################################

## whole-ROI annotated maps
p_bt288 <- ImageDimPlot(xenium, fov = "ROI_BT288_1", group.by = "annotation_v2",
                        cols = parade_map, size = 1, border.size = NA,
                        dark.background = TRUE, nmols = 20000) + NoAxes()
p_tma8  <- ImageDimPlot(TMA8,  fov = "ROI_TMA8",  group.by = "annotation_v2",
                        cols = parade_map, size = 1, border.size = NA,
                        dark.background = TRUE, nmols = 20000) + NoAxes()
p_tma11 <- ImageDimPlot(TMA11, fov = "ROI_TMA11", group.by = "annotation_v2",
                        cols = parade_map, size = 1, border.size = NA,
                        dark.background = TRUE, nmols = 20000) + NoAxes()

## export one PDF per ROI
outdir <- "All_ROI_ImageDimPlots"; dir.create(outdir, showWarnings = FALSE)
for (nm in names(xenium@images)) {
  p <- ImageDimPlot(xenium, fov = nm, group.by = "annotation_v2", cols = parade_map,
                    size = 0.5, border.size = NA, dark.background = TRUE,
                    nmols = 20000) + NoAxes() + NoLegend()
  ggsave(file.path(outdir, paste0(nm, ".pdf")), p, width = 6, height = 6, bg = "black")
}

################################################################################
## 7. Marker feature maps (ameloblast / junctional / sulcus)
################################################################################

pu_or <- rev(brewer.pal(n = 8, name = "PuOr"))

## ameloblast enamel markers
ImageFeaturePlot(TMA8,  fov = "ROI_TMA8",  features = c("KLK4","MMP20","ENAM","AMBN"),
                 size = 0.75, border.size = NA, cols = pu_or)
ImageFeaturePlot(xenium, fov = "ROI_BT288_1", features = c("ALPL","KRT6C","KRT6B","AMBN"),
                 size = 0.75, border.size = NA, cols = pu_or)
## junctional / sulcus markers
ImageFeaturePlot(TMA11, fov = "ROI_TMA11", features = c("SPRR1B","SPRR2F","AMTN","AMBN"),
                 size = 0.75, border.size = NA, cols = pu_or)

################################################################################
## 8. General oral / dental module scores (spatial)
################################################################################

General_oral   <- c("CRYAB","NPPC","ITM2A","PCBD1","ZBED2","TUSC1","FXYD3","GCSH",
                    "DCTN3","TNFSF10","IFITM3","SELENOP","IRX1","CXCL14","DUSP1")
General_dental <- c("DISC1","AUTS2","ADGRB3","FBN2","TSHZ2","EYA1","PDE4D","EPHA7",
                    "ADGRV1","MTUS2","CUX2","MEG8","GRAMD1B","VWDE","FMN1","DCC",
                    "PCAT1","RIMS2","ENOX1")
map_genes <- function(g) { rn <- rownames(xenium[[DefaultAssay(xenium)]])
                           unique(rn[match(toupper(g), toupper(rn))][!is.na(rn[match(toupper(g), toupper(rn))])]) }

xenium <- AddModuleScore(xenium, features = list(map_genes(General_oral)),
                         name = "General_oral_score",  ctrl = 5, nbin = 12)
xenium <- AddModuleScore(xenium, features = list(map_genes(General_dental)),
                         name = "General_dental_score", ctrl = 5, nbin = 12)

ImageFeaturePlot(xenium, fov = "ROI_BT288_1", features = "General_dental_score1",
                 size = 0.75, border.size = NA, cols = pu_or)
ImageFeaturePlot(xenium, fov = "ROI_BT288_1", features = "General_oral_score1",
                 size = 0.75, border.size = NA, cols = pu_or)

################################################################################
## 9. Supplementary marker dot plot
################################################################################

genelist <- c("SPINK5","EHF","PAPPA","DLX5","FBN2","SCUBE3","LGR6","LEF1","CHGB","AXIN2",
              "NOTUM","SOX2","KRT6C","KRT6B","ALPL","SP6","VWDE","KIF5C","ENAM","AMBN",
              "GJB6","GJB2","ODAM","AMTN","SPRR1B","SPRR2F","CCL2","COL3A1","MKI67",
              "PTPRC","S100B","VWF")
genes_present <- genelist[genelist %in% rownames(xenium)]

p <- DotPlot_scCustom(xenium, features = genes_present, group.by = "annotation_v2")
ggplot(p$data, aes(y = id, x = features.plot, size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, colour = "black") +
  theme_bw() + scale_size_area(max_size = 10) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "red3", limits = c(-2.5, 2.5)) +
  labs(fill = "Relative expression", size = "Percent expressed", x = "Cluster", y = "Gene") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position = "bottom")

################################################################################
## NOTE
## - Immune subset re-clustering + annotation_immune: see immune subset script.
## - Spatial neighbourhood enrichment (Squidpy-equivalent, developmental pairs and
##   immune vs oral/dental): see spatial_adjacency.R.
################################################################################



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

