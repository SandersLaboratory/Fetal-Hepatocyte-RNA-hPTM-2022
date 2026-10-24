# PCA — README

Principal component analysis plots for the Adult–Fetal project (RNA-seq and hPTM data), supporting both 2D (`ggfortify::autoplot`) and 3D (`plot3D::scatter3D`) styles. All four full scripts share the same template: user settings at the top, edgeR TMM normalization (RNA-seq) or direct input (hPTM), `prcomp(..., scale.=TRUE)`, then plotting.

## `PCA.R`
- **Purpose:** PCA of LCM RNA-seq counts, colored by Adult/Fetal.
- **Inputs:** `AF_Genewiz_all_counts_no246.csv` + `afg_metadata_no246.csv` in `Adult-Fetal/Data Files/rna-seq/AF_Genewiz/` (lines 15–18).
- **Outputs:** Plot prints to the viewer; no files written (a commented-out block would save TMM-normalized CPMs — leave disabled, per the `#DO NOT USE heatmaply(pcadata)` note).
- **Dependencies:** `tidyverse`, `edgeR`, `plot3D`, `ggfortify`, `directlabels`, `rstudioapi`.
- **Usage/Notes:** Settings at lines 14–66: `samplesToGraph = "Metadata"`, `metadataColor = "adultfetal"` (Adult green4 / Fetal red), `metadataShape` currently `"adultfetal"` (both triangle; a replicate-shape block is overridden below it), `typeOfPCA = "3D"`, view angles `leftToRight = -30`, `bottomToTop = 40`, `pointSize = 2.5`. Switch `typeOfPCA` to `"2D"` and set `styleOfPCA` to `"Shape"`/`"Label"`, `PCX`/`PCY`, `lines = TRUE` for replicate-connecting lines.

## `PCA_hPTM.R`
- **Purpose:** PCA over merged hPTM mass-spec tables (cell + colony data; liver optional), then Euclidean distances between group centroids in PC space.
- **Inputs:** `hPTM_Cells - 020716 - 2022.csv`, `Colonies - 012320.csv`, optionally `Liver - 080919.csv`, and `hptm_metadata_2022.csv` in `Adult-Fetal/Data Files/proteomics/hPTM/` (lines ~68–90).
- **Outputs:** 3D (or 2D) plot to the viewer; two printed distance matrices: (1) 4×4 between Adult/Fetal × Cells/HostCol centroids, (2) 6×6 between subgroup centroids (A-20, F-20, FEHS, FEHD, ADH, FEHL) — computed on PCs 1–3.
- **Dependencies:** `tidyverse`, `edgeR`, `plot3D`, `ggfortify`, `directlabels`, `rstudioapi`.
- **Usage/Notes:** `includeLiver = FALSE` (line ~58) currently drops the 6 liver samples and their metadata rows; set `TRUE` to include them. Color = adultfetal, shape = `Group` (A-20/F-20/FEHS/FEHD/ADH/FEHL). Note the header comment still says "RNA-Seq Data" — this copy was adapted for hPTM.

## `PCA_cells.R`
- **Purpose:** PCA of AF cell RNA-seq, LCM RNA-seq, or both combined — with ComBat-seq batch correction when combining.
- **Inputs:** `AF_Cells_raw_counts.csv`, `AF_Genewiz_all_counts_no246.csv`, `AF_Merged_raw_counts.csv` and matching metadata files (`AF_Cells_metadata.csv`, `afg_metadata_no246.csv`, `AF_Merged_metadata.csv`) in `Adult-Fetal/Data Files/rna-seq/AF_Genewiz/` (lines ~102–126).
- **Outputs:** Plot to the viewer; no files written (the merge/normalized-count `write.csv` lines are commented out).
- **Dependencies:** `tidyverse`, `edgeR`, `plot3D`, `ggfortify`, `directlabels`, `sva` (only when batch-correcting), `org.Hs.eg.db`, `org.Rn.eg.db`, `annotate`, `rstudioapi`.
- **Usage/Notes:**
  - `inputToGraph = "Both"` (line 16) selects Cells / LCM / Both. When `"Both"`, `comBatSeq = TRUE` and `ComBat_seq` corrects on the metadata `batch` column (no biological covariate protected — see the SVA README note).
  - A block that rebuilds `AF_Merged_raw_counts.csv` from the two count tables is commented out (lines ~107–114); the script currently reads the pre-built merged file.
  - The 2D "Shape" branch applies a legend-free figure style (no legend, no axis titles, large tick labels) — tuned for posters.
  - Ends with John's top-PC-contributor snippet (same code as `John Code - Top PCA Contributors.R`).

## `John Code - Top PCA Contributors.R`
- **Purpose:** Rank genes by their contribution to one principal component.
- **Inputs:** None directly — expects `cpmdata` (genes × samples CPM matrix) and `groups` (sample metadata) to already exist in the R session (e.g. from running one of the PCA scripts above; `cpmdata` is only a comment — `PCA_cells.R` runs the same lines on `pcadata`).
- **Outputs:** `percents` — a named vector of (loading²) per-gene PC2 variance contributions, sorted descending.
- **Dependencies:** None beyond base R (`prcomp` from `stats`).
- **Usage/Notes:** Not standalone — source it after a PCA workflow, or copy the lines to the end of a script. Change `contribution[,"PC2"]` to target a different PC.

## Duplicates of GenericCode versions
- `CalculatePCACentroidEuclidianDistance.R` — identical to `GenericCode/R Scripts/PCA/CalculatePCACentroidEuclidianDistance.R` — see that README.
- `PCA_Proteomics.R` — identical to `GenericCode/R Scripts/PCA/PCA_Proteomics.R` — see that README.
