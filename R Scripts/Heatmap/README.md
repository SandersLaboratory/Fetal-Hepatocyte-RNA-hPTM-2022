# Heatmap — README

Interactive `heatmaply` figures of histone-PTM differential results for the Adult–Fetal project.

## `Heatmap Histone PTMs.R`
- **Purpose:** Build three interactive heatmaps from the all-study hPTM comparison table, plus (leftover) grouped bar-chart code for the 13 significant marks.
- **Inputs:**
  - `hPTM All Study Comparison - Merged 2022.csv` in `Adult-Fetal/Data Files/BarGraph/` (line 5; rows = hPTM marks, columns = effect + q-value pairs for Colonies, Liver, LAP+/OC2+, LAP+/OC2−, LAP+).
  - `hPTM - Significant Table7_Feb2021.csv` (same directory) for the bar-chart section.
- **Outputs (working directory):**
  - `hPTM_heatmaply_cor.pdf` — heatmap of only significant marks, point size = −log10(q), "ns" cell notes for non-significant cells.
  - `hPTM_heatmaply_cor_full.png` — all marks, with `*` (q ≤ 0.05) / `^` (0.05 < q ≤ 0.1) cell annotations.
  - `hPTM_heatmaply.png` — plain heatmap of significant marks, blue–white–red scale.
- **Dependencies:** `tidyverse`, `data.table`, `RColorBrewer`, `ggpattern`, `plotly`, `heatmaply` (lines 23–28), plus `rstudioapi` for working-directory auto-detection.
- **Usage/Notes:**
  - Set `significance` (line 7, currently 0.05) and `dropLiver` (line 10; when `TRUE`, the Liver columns are dropped from the heatmaps — currently `FALSE`).
  - The bottom ~150 lines are bar-chart code duplicated from the BarGraph scripts (mean ± SE of the 13 significant marks) — not wired to the heatmap workflow; safe to ignore or delete.
  - Contains a workaround for the `hclust` NA/NaN error via `complete.cases` filtering before clustering.

## `Heatmap&Dotplot_hPTM-logFC.R`
- **Purpose:** Heatmap + dot plot of logFC values for the overlapping significant hPTMs across the three fetal-vs-adult comparisons: Fetal Colony, Dual-Positive (LAP+/OC2+), Single-Positive (LAP+/OC2−).
- **Inputs:** `hPTM All Study Comparison -Final Overlap.csv` in `Adult-Fetal/Data Files/proteomics/` (line 18; rows = hPTM marks, alternating logFC / q-value columns).
- **Outputs:** Plots print to the RStudio viewer only — no files are written (the commented-out `ggsave` line is absent; `heatmaply` and `plotOverlap` render interactively).
- **Dependencies:** `tidyverse`, `plotly`, `heatmaply`, `data.table`, `org.Hs.eg.db` (loaded but effectively unused in this script).
- **Usage/Notes:**
  - The script first unloads all previously loaded packages (lines 23–28) — source it in a fresh session or it will detach your other packages.
  - Heatmap: logFC columns 1/3/5, blue–white–red diverging scale clamped to ±30.
  - Dot plot (`plotOverlap`): y = comparison, x = hPTM (flipped coords); fill = logFC (blue–white–red, ±30), point size = q value (reverse-scaled so smaller q = bigger dot), black-stroked circles.
  - Authored by Anders Ohman, Sanders Lab, Spring 2020 (header comment).
