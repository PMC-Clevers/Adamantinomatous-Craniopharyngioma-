################################################################################
## cellchat_primary10x.R
## Ligand-receptor / cell-cell communication analysis (CellChat) on the primary
## ACP scRNA-seq cohort, immune + epithelial fine cell types.
## CellChat v2.2.0.9001.
################################################################################

library(CellChat)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
options(stringsAsFactors = FALSE)

## =============================================================================
## 1. Build the CellChat object from the annotated Seurat object
##    (fine annotation: immune subtypes + epithelial sub cell types)
## =============================================================================

all_prim_cc <- readRDS("all_prim_with_fine_annotation.rds")
Idents(all_prim_cc) <- "finer_annotation"

DefaultAssay(all_prim_cc) <- "RNA"
data.input <- GetAssayData(all_prim_cc, assay = "RNA", layer = "data")   # log-normalized
meta <- data.frame(labels = Idents(all_prim_cc), row.names = colnames(all_prim_cc))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")

## interaction database (human): secreted signalling, ECM-receptor, cell-cell contact
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

## =============================================================================
## 2. Core communication inference
## =============================================================================

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, population.size = TRUE)   # weight by group size
cellchat <- filterCommunication(cellchat, min.cells = 10)        # >= 10 cells per group
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, "cellchat_ACP_primary10x.rds")

## =============================================================================
## 3. Full ligand-receptor table (df.net) - export for filtering
## =============================================================================

df.net <- subsetCommunication(cellchat)          # all significant L-R pairs
write.csv(df.net, "CellChat_ACP_LR_network_full.csv", row.names = FALSE)

## =============================================================================
## 4. Global network (Supplementary Figure S11)
## =============================================================================

## A: total aggregated interactions (chord)
pdf("cellchat_TOTAL_network.pdf", width = 12, height = 11); par(xpd = TRUE)
netVisual_aggregate(cellchat, signaling = NULL, layout = "chord",
                    vertex.weight = as.numeric(table(cellchat@idents)))
dev.off()

## B: outgoing / incoming signalling-role heatmaps
ph_out <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 12, height = 18)
ph_in  <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 12, height = 18)
pdf("cellchat_signalingRole_heatmap.pdf", width = 20, height = 18)
draw(ph_out + ph_in)
dev.off()

## C: representative developmental pathway chords
for (pw in c("FGF","BMP")) {
  if (pw %in% cellchat@netP$pathways) {
    pdf(paste0("cellchat_", pw, "_chord.pdf"), width = 9, height = 9); par(xpd = TRUE)
    netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
    dev.off()
  }
}

## =============================================================================
## 5. GAS6-AXL/MERTK: oral/gingival -> myeloid (Supplementary Figure S10D-E)
## =============================================================================

gingival <- c("Basal progenitors","CCL2+ Keratinocyte","Junction Basal/Keratinocyte",
              "Pre-Mesenchymal-like","Stress","B/Plasma cells")
myeloid  <- c("Macrophage M1","Macrophage M2","TAM","cDC2")

## D: GAS chord (senders -> myeloid receivers)
if ("GAS" %in% cellchat@netP$pathways) {
  pdf("cellchat_GAS_chord.pdf", width = 8, height = 8); par(xpd = TRUE)
  netVisual_chord_gene(cellchat, sources.use = gingival, targets.use = myeloid,
                       signaling = "GAS", lab.cex = 0.7, legend.pos.x = 8, legend.pos.y = 30)
  dev.off()
}

## E: GAS6 ligand + TAM-receptor expression dot plot
DefaultAssay(all_prim_cc) <- "RNA"
DotPlot(all_prim_cc, features = c("GAS6","PROS1","AXL","MERTK","TYRO3"),
        group.by = "finer_annotation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("GAS6_TAMreceptor_dotplot.pdf", width = 7, height = 10)

## GAS interactions as a table
gas <- df.net[df.net$pathway_name == "GAS", ]
gas <- gas[order(-gas$prob), ]
write.csv(gas, "GAS6_AXL_interactions.csv", row.names = FALSE)

## =============================================================================
## 6. Immune-directed signalling bubbles (epithelial -> immune)
## =============================================================================

immune_targets <- c("Macrophage M1","Macrophage M2","TAM","cDC2",
                    "CD4+ T","CD8+ T","Tregs","B/Plasma cells","Monocytes/Neutrophils")
epi_sources    <- setdiff(levels(cellchat@idents), immune_targets)

pdf("cellchat_epithelial_to_immune_bubble.pdf", width = 12, height = 16)
netVisual_bubble(cellchat, sources.use = epi_sources, targets.use = immune_targets,
                 remove.isolate = TRUE, angle.x = 45)
dev.off()
