################################################################################
## 01_organoid_pipeline_acp.R
## ACP patient-derived organoid scRNA-seq pipeline.
## Organoid subset -> QC -> SCT -> integration (CCA + RPCA) -> clustering ->
## cell-cycle scoring -> annotation (annotation_orgs_v1) -> marker plots.
## RPCA is the integration used for the final clustering/annotation.
################################################################################

library(Seurat)
library(SeuratWrappers)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(plyr)          # mapvalues
library(RColorBrewer)
library(future)
options(future.globals.maxSize = 60 * 1024^3)

## =============================================================================
## 1. Subset organoids from the full object and QC
## =============================================================================

acp <- readRDS("20250424_acp_unprocessed_all.rds")

all_orgs <- subset(acp, subset = test_genotype %in%
  c("BT322_Org","BT506_Org","BT259_Org","BT325_Org","BT154_Org","BT288_Org","BT155_Org"))
table(all_orgs$test_genotype)

all_orgs <- subset(all_orgs, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

## =============================================================================
## 2. Normalize, PCA, and unintegrated embedding (inspection)
## =============================================================================

all_orgs[["RNA"]] <- JoinLayers(all_orgs[["RNA"]])
all_orgs[["RNA"]] <- split(all_orgs[["RNA"]], f = all_orgs$test_genotype)
all_orgs <- SCTransform(all_orgs)
all_orgs <- RunPCA(all_orgs)
ElbowPlot(all_orgs, ndims = 40)

all_orgs <- FindNeighbors(all_orgs, reduction = "pca", dims = 1:25, nn.eps = 0.5)
all_orgs <- FindClusters(all_orgs, resolution = 0.5, n.start = 10)
all_orgs <- RunUMAP(all_orgs, dims = 1:25, min.dist = 0.5, reduction.name = "umap.unintegrated")
DimPlot(all_orgs, label = TRUE, split.by = "test_genotype", reduction = "umap.unintegrated")

## =============================================================================
## 3. CCA integration (comparison)
## =============================================================================

all_orgs <- IntegrateLayers(all_orgs, method = CCAIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.cca",
                            normalization.method = "SCT", verbose = FALSE)
all_orgs[["RNA"]] <- JoinLayers(all_orgs[["RNA"]])
all_orgs <- FindNeighbors(all_orgs, reduction = "integrated.cca", dims = 1:20)
all_orgs <- FindClusters(all_orgs, resolution = 0.4, cluster.name = "cca_clusters")
all_orgs <- RunUMAP(all_orgs, reduction = "integrated.cca", dims = 1:20,
                    min.dist = 0.1, reduction.name = "umap.ccaintegrated")
DimPlot(all_orgs, label = TRUE, reduction = "umap.ccaintegrated", split.by = "test_genotype")

## =============================================================================
## 4. RPCA integration (final: used for clustering + annotation)
## =============================================================================

all_orgs[["RNA"]] <- split(all_orgs[["RNA"]], f = all_orgs$test_genotype)
all_orgs <- IntegrateLayers(all_orgs, method = RPCAIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.rpca", dims = 1:30,
                            normalization.method = "SCT", verbose = TRUE)
all_orgs[["RNA"]] <- JoinLayers(all_orgs[["RNA"]])

all_orgs <- FindNeighbors(all_orgs, reduction = "integrated.rpca", dims = 1:20)
all_orgs <- FindClusters(all_orgs, resolution = 0.5, cluster.name = "rpca_clusters20")
all_orgs <- RunUMAP(all_orgs, reduction = "integrated.rpca", dims = 1:20,
                    min.dist = 0.1, reduction.name = "umap.rpca", return.model = TRUE)
DimPlot(all_orgs, label = TRUE, split.by = "test_genotype", reduction = "umap.rpca")

## =============================================================================
## 5. Markers and cell-cycle scoring
## =============================================================================

all_orgs <- PrepSCTFindMarkers(all_orgs, assay = "SCT", verbose = TRUE)
all_orgs.markers <- FindAllMarkers(all_orgs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all_orgs.markers, "allorgs_markers.csv")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all_orgs <- CellCycleScoring(all_orgs, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(all_orgs, features = c("PCNA","TOP2A","MCM6","MKI67"), ncol = 2)
FeaturePlot(all_orgs, features = c("S.Score","G2M.Score"), reduction = "umap.rpca") &
  scale_color_gradientn(colors = rev(brewer.pal(n = 8, name = "PuOr"))) & NoAxes()

## =============================================================================
## 6. Annotation (annotation_orgs_v1)
## =============================================================================

Idents(all_orgs) <- "rpca_clusters20"
cluster.ids_v1 <- c(
  "CCK/DKK1/FST",          # 0
  "LGR6+ CL-like",         # 1
  "Cluster Cells",         # 2
  "CCL2/IRX1 OE",          # 3
  "Cycling",               # 4
  "Intermediate",          # 5
  "Gingival epithelium",   # 6
  "OEE/SR/SI",             # 7
  "Cycling",               # 8
  "CCL2+ keratinocytes",   # 9
  "HERS/ERM-like",         # 10
  "Pre-cluster",           # 11
  "Cycling",               # 12
  "Pre-ameloblast/IEE",    # 13
  "PDL-like"               # 14
)
names(cluster.ids_v1) <- as.character(0:14)
all_orgs$annotation_orgs_v1 <- plyr::mapvalues(as.character(Idents(all_orgs)),
                                               from = names(cluster.ids_v1), to = cluster.ids_v1)
Idents(all_orgs) <- "annotation_orgs_v1"

DimPlot(all_orgs, reduction = "umap.rpca", group.by = "annotation_orgs_v1",
        label = TRUE, label.size = 6) + NoLegend()
saveRDS(all_orgs, "all_orgs_rpca_annotated.rds")

## =============================================================================
## 7. Marker dot plot (annotation supplementary)
## =============================================================================

markerlist <- c(
  "COL3A1","TAGLN",                       # PDL-like
  "SPINK5","HOPX","ODAM","GJB6","ATF3",   # gingival epithelium
  "DKK4","HUNK","AXIN2","SP6","AMBN",     # cluster cells / ameloblast
  "HMGA2","LAMC1",                        # HERS / ERM
  "FBN2","AFF3","RIMS2","SOX6","SOX5","CLDN4",  # OEE / SI
  "KIF5C","UNC5C","ETV5",                 # IEE
  "LGR6","FRZB",                          # CL-like
  "CCL2","CXCL1",                         # CCL2 module
  "MKI67")                                # cycling

DotPlot_scCustom(all_orgs, features = markerlist, group.by = "annotation_orgs_v1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ameloblast / lineage feature plots
FeaturePlot_scCustom(all_orgs, features = c("AMBN","SP6","VWDE","KIF5C"),
                     num_columns = 2, reduction = "umap.rpca") & theme_void()

## =============================================================================
## 8. General oral / dental module scores
## =============================================================================

General_oral <- c("CRYAB","NPPC","ITM2A","PCBD1","ZBED2","TUSC1","FXYD3","GCSH",
                  "DCTN3","TNFSF10","IFITM3","SELENOP","IRX1","CXCL14","DUSP1")
General_dental <- c("DISC1","AUTS2","ADGRB3","TSHZ2","EYA1","SDK1","PDE4D","EPHA7",
                    "ADGRV1","MTUS2","CUX2","MEG8","SEMA6D","RBFOX1","KMT2C","TEAD1",
                    "YAP1","GRAMD1B","VWDE","FMN1","DCC","PCAT1","RIMS2","ENOX1",
                    "PCDH15","PTPRT","AGAP1")

add_two_panels <- function(obj, oral_genes, dental_genes, ctrl = 5, nbin = 12) {
  rn <- rownames(obj[[DefaultAssay(obj)]])
  map_genes <- function(g) { h <- rn[match(toupper(g), toupper(rn))]; unique(h[!is.na(h)]) }
  obj <- AddModuleScore(obj, features = list(map_genes(oral_genes)),
                        name = "General_oral_score", ctrl = ctrl, nbin = nbin)
  obj$General_oral <- obj$General_oral_score1
  obj <- AddModuleScore(obj, features = list(map_genes(dental_genes)),
                        name = "General_dental_score", ctrl = ctrl, nbin = nbin)
  obj$General_dental <- obj$General_dental_score1
  obj
}
all_orgs <- add_two_panels(all_orgs, General_oral, General_dental)

FeaturePlot(all_orgs, features = c("General_oral_score1","General_dental_score1"),
            order = TRUE, reduction = "umap.rpca") &
  scale_color_gradientn(colors = rev(brewer.pal(n = 8, name = "PuOr"))) & NoAxes()
