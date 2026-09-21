# Adamantinomatous Craniopharyngioma: developmental single-cell and spatial analysis

Code accompanying *"Adamantinomatous Craniopharyngioma, a benign pituitary tumor
recapitulating the epithelial program of tooth and gingival development."*

This repository contains the analysis scripts for the single-cell RNA-seq
(primary tumor and organoids), the Xenium spatial transcriptomics, the immune
subtyping, the CellChat ligand-receptor analysis, and the Squidpy-equivalent
spatial neighborhood-enrichment analysis.

 - **Primary tumor scRNA-seq** (7 patients) to resolve the epithelial cell states
  and map them onto dental and oral/gingival developmental references.
- **Patient-derived organoids** (7 lines), including directed differentiation
  toward ameloblast- and periodontal-ligament-like fates, to test whether
  CTNNB1-mutant ACP cells retain this developmental potential.
- **Xenium spatial transcriptomics** (5 tumors, 15 ROIs) to test whether the
  cell states are spatially organized along the developmental axis.
  
## Environment

- R 4.6.1
- Seurat 5.5.1, SeuratWrappers, scCustomize
- harmony, symphony 0.1.3, monocle3, CellChat 2.2.0.9001
- FNN, pheatmap, ggplot2, dplyr, purrr
- SoupOrCell (genotype demultiplexing, run upstream)

Processed Seurat objects and DESeq2 outputs are deposited on Zenodo (DOI in the
manuscript). Analysis code is version-tagged in this repository.

## Scripts

### `00_full_pipeline_acp.R` — primary 10x scRNA-seq
Per-lane 10x load and QC, SoupOrCell SNP genotype demultiplexing (inter-patient
doublet / unassigned removal), merge, and RPCA integration of the full object.
The primary-tissue subset is re-integrated (`integrated.rpca100`), and
non-epithelial cells are removed in two passes: (1) immune + endothelial
clusters, re-integrate, then (2) the glial cluster that separates on
re-clustering, re-integrate again. Both `FindNeighbors` and `RunUMAP` are run on
the same fresh `integrated.rpca` so clustering and embedding share one space.
The clean epithelial object is split into dental and oral/gingival compartments.
Intra-genotype doublets were not separately assessed (stated as a limitation).

### `01_organoid_pipeline_acp.R` — organoid scRNA-seq
Organoid subset, QC, SCTransform, CCA and RPCA integration (RPCA used for the
final clustering), cell-cycle scoring, annotation (`annotation_orgs_v1`), and
marker / oral-dental module-score plots.

### Reference mapping (out-group controls, Reviewer 3.1)
Projection of the ACP dental and oral epithelial subsets onto dental, gingival,
epidermis and esophagus references using two frameworks: Seurat
`FindTransferAnchors`/`MapQuery` (with `MappingScore`), and Symphony
(`build_symphony_ref_manual` → `mapQuery` → `knnPredict`, k = 5) on a joint
four-tissue reference.

### Pseudotime (Monocle3)
Pseudotime on the dental-like compartment, rooted in the Transitioning
population; `graph_test` for stage markers and `plot_genes_in_pseudotime` along
the Transitioning → Ameloblast axis.

### `xenium_analysis_final.R` — Xenium spatial transcriptomics
Load runs (TMA + BT288, BT259 excluded), QC (`nCount_Xenium > 10`), Harmony
integration by run, clustering, `annotation_v2`, ROI splitting (13 TMA cores +
2 BT288 pieces = 15 ROIs), and the paper figures (spatial cluster maps, marker
feature plots, oral/dental module scores, supplementary dot plot).

### `xenium_immune_subset.R` — Xenium immune subtyping
Subclustering and annotation of the myeloid/immune compartment into six subtypes
(TAM GPNMB+/APOE+, Macrophages, Macrophages CCL2+, Macrophages CD86+, T-cells
CD8, T-cells CD4/Treg), merged back into the full object as `annotation_immune`
(`xenium_with_immune_subtypes.rds`).

### CellChat — ligand-receptor analysis (primary 10x)
Cell-cell communication on the primary scRNA-seq cohort (human CellChatDB,
`computeCommunProb` with `population.size = TRUE`, `subsetCommunication`),
including the GAS6-AXL/MERTK oral/gingival → myeloid axis and the shared
MIF / midkine recruitment signals.

### `spatial_adjacency.R` — Squidpy-equivalent spatial neighborhood enrichment
Neighborhood-enrichment z-scores from an undirected, symmetrised k-nearest-
neighbour graph (k = 6) with 1,000 label permutations,
`z = (observed - mean(permuted)) / SD(permuted)`, computed within each ROI so no
edges cross tissue pieces. Includes:
1. developmental cell-type pairs (adjacent vs non-adjacent controls),
2. developmental stages (sequential vs stage-skipping),
3. the unbiased z-score heatmap across all populations,
4. the immune compartment: each immune subtype vs the collapsed oral and dental
   fractions (`annotation_immune`).
Cohort-level significance per pair uses a one-sample Wilcoxon test with
Benjamini-Hochberg correction. Z-scores are aggregated per patient (median
across ROIs); pairs present in fewer than three patients are excluded.

**Positive controls** (adjacent, expected to be spatially close):
`OE → TR → DE → DEOEE → SI/SR / Cluster Cells → pre-AMB → AMB → JK/JB`

**Negative controls** (non-adjacent, expected to be spatially distant):
`OE → SI/SR (skip cap)`, `OE → AMB (skip all)`, `DE → AMB (skip)`, `TR → AMB (skip)`

## Methods summary

### Nearest-neighbour distance analysis
For each pair of cell populations, the median nearest-neighbour distance (µm) from all cells of type A to the closest cell of type B was computed using the `FNN` package (k = 1). Statistical significance was assessed by 1,000-iteration label permutation testing: cell type labels were shuffled within each ROI while preserving spatial coordinates, generating a null distribution of expected distances. The p-value was defined as the fraction of permuted distances ≥ the observed distance (low p-value = cells closer than expected by chance). Permutations were parallelized across available CPU cores.

### Neighborhood enrichment z-score
An R implementation equivalent to the Squidpy `gr.nhood_enrichment` function (Palla et al., *Nature Methods* 2022) was used. For each patient (ROIs pooled), a k = 6 spatial nearest-neighbor graph was built using physical cell coordinates. Co-occurrence counts between cell type pairs were computed on the undirected (symmetrized) graph and compared against 1,000 label permutations, yielding z-scores: z = (observed − mean(null)) / sd(null). Positive z-scores indicate enriched spatial co-occurrence; negative z-scores indicate depletion.

### Cell type reliability filter
Only cell types present in ≥ 3 of 5 patients (with ≥ 20 cells per patient) were included in the unbiased proximity landscape analyses.

## Requirements

```r
# R >= 4.2
install.packages(c(
  "Seurat",       # >= 5.0
  "dplyr",
  "ggplot2",
  "parallel",
  "purrr",
  "tidyr",
  "tibble",
  "FNN",
  "pheatmap",
  "viridis"
))
```

## Input data

| File | Description |
|------|-------------|
| `20251103_xenium-nobt259.rds` | Seurat object with Xenium spatial transcriptomics data (5 ACP patients, 14 ROIs) |

The Seurat object requires:
- `@meta.data$annotation_v2` — cell type annotation column
- `@images` — named FOV objects per ROI with spatial coordinates

## Usage

```r
# Set working directory to folder containing the xenium RDS file
setwd("path/to/data")

# Run full analysis (takes ~2-4 hours depending on CPU cores)
source("spatial_adjacency_v4.R")
```

Results are saved to `results/spatial_adjacency_v4/`.

**Note:** Neighborhood enrichment z-scores are cached per patient as `zscore_k6_ACP00X.rds`. If these files exist, they are loaded rather than recomputed, making reruns fast.


## Patient cohort

| Patient | Xenium profiling | ROIs | Notes |
|---------|-----------------|------|-------|
| ACP004  | Whole slide (BT288) | 1 | 401,841 cells |
| ACP005  | TMA | 2 | |
| ACP006  | TMA | 3 | |
| ACP007  | TMA | 7 | |
| ACP008  | TMA | 1 | |

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_CELLS` | 20 | Minimum cells per population per ROI |
| `N_PERM` | 1,000 | Permutation iterations |
| `SEED` | 42 | Random seed |
| `K_NHOOD` | 6 | k for neighborhood enrichment (Squidpy default) |
| Min patients | 3 | Minimum patients for proximity landscape |

## Citation

If you use this code, please cite:

> Verweij L, et al. (2026). Adamantinomatous Craniopharyngioma, a benign pituitary tumor recapitulating development of tooth and gingiva. *Neuro-Oncology*. doi: [XXX]

And the key dependencies:

> Palla G, et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178.

> Hao Y, et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42, 293–304.

## License

MIT License — see `LICENSE` for details.

## Data availability
Processed Seurat objects, annotations and DESeq2 differential-expression outputs
are deposited on Zenodo (versioned DOI, see manuscript). Xenium run-level quality
metrics are provided in Supplementary Table S6.

## Contact

Laurens Verweij — Prinses Máxima Centrum / Utrecht University
Supervisors: Hans Clevers, Marc van de Wetering
'

