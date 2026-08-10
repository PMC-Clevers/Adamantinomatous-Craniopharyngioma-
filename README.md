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

## Output files

### CSV results
| File | Description |
|------|-------------|
| `spatial_nn_individual_per_ROI.csv` | Permutation test results per ROI per pair |
| `spatial_nn_individual_per_patient.csv` | Aggregated per patient |
| `spatial_nn_stages_per_patient.csv` | Stage-level results (pooled ROIs) |
| `all_pairwise_distances.csv` | Raw pairwise NN distances — all cell types × all patients |
| `pairwise_distance_summary.csv` | Median pairwise distances across patients (≥3 patients filter) |
| `full_pairwise_zscore_k6_squidpy.csv` | Full z-score matrix per patient |
| `zscore_summary_k6_squidpy.csv` | Median z-scores across patients (≥3 patients filter) |

### Key figures
| File | Description |
|------|-------------|
| `plot1_individual_violin_v2` | NN distance distributions — adjacent vs non-adjacent pairs |
| `plot2_individual_pvalue_heatmap` | Permutation test p-values per patient × pair |
| `plot3_stage_pvalue_heatmap` | Stage-level permutation test p-values |
| `plot3_stage_violin` | Stage-level distance distributions |
| `plot4_stage_progression` | Per-patient Bud → Cap → Bell → Eruption progression |
| `plot5_distance_heatmap_capped_purple_orange` | Hierarchical clustering — NN distance |
| `plot6_zscore_heatmap_hierarchical` | Hierarchical clustering — z-score |
| `plot_density_per_patient_updated` | Distance density distributions per patient |
| `plot_ranked_zscore_updated` | Ranked neighborhood enrichment z-scores |
| `distance_landscape_combined` | Unbiased proximity landscape — all cell types (distance) |
| `zscore_landscape_full` | Unbiased proximity landscape — all cell types (z-score) |
| `distance_[CellType].pdf` | Per cell type distance landscape (25 files) |
| `zscore_k6_[CellType].pdf` | Per cell type z-score landscape (25 files) |

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

## Contact

Laurens Verweij — Prinses Máxima Centrum / Utrecht University
Supervisors: Hans Clevers, Marc van de Wetering
'

