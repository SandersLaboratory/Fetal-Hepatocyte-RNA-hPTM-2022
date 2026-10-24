# Volcano — README

Differential-expression (edgeR) + volcano-plot pipelines for the Adult–Fetal RNA-seq data. All scripts follow the same skeleton — User Settings → edgeR GLM (`estimateGLM*Disp` → `glmFit` → `glmLRT`) → `EnhancedVolcano` plot → CSV/PNG outputs — with different contrasts, datasets, and annotation styles. Run from RStudio: every script derives the repo root from its own path via `rstudioapi` (e.g. `Volcano_2.0.R` lines 51–59).

## `Volcano_2.0.R`
- **Purpose:** Full DE pipeline + publication-style volcano for the LCM dataset (Fetal vs Adult), with three selectable plot styles.
- **Inputs:** `AF_Genewiz_all_counts_no246.csv` + `afg_metadata_no246.csv` in `Adult-Fetal/Data Files/rna-seq/AF_Genewiz/` (lines 15–16); `comparisonGroup = "adultfetal"`, `inputContrast = "Fetal-Adult"` (lines 18–19).
- **Outputs** (input directory, lines 263–274): `DEAnalysis_adultfetal_Fetal-Adult.csv` (all genes + geneid/genesymbol), `SigUp_…csv` / `SigDown_…csv` (FDR < 0.05), and `<style>_adultfetal_Fetal-Adult.png`.
- **Dependencies:** `tidyverse`, `edgeR`, `annotate`, `EnhancedVolcano`, `org.Hs.eg.db`, `org.Rn.eg.db`.
- **Usage/Notes:**
  - `style = "Minimal"` (current; figure-ready: no titles/labels, red = up, royalblue = down, black = ns, fixed axes ±6 / 0–6); `"Simple"` (red = significant); `"Enhanced"` (labeled points, FC cutoff 1.0). `userFDR = 0.05`.
  - Gene symbols are mapped ENSEMBL → Entrez → SYMBOL via `org.Rn.eg` lookups; unmapped rows get a blank symbol.
  - A commented-out tail (lines ~276+) would prune the count matrix to significant genes for a filtered-counts export.

## `Volcano_Cells.R`
- **Purpose:** Same DE + volcano pipeline as `Volcano_2.0.R`, but switchable between three datasets.
- **Inputs:** Set `comparison = "LCM"` / `"Dual"` / `"Single"` (line 18); picks LCM counts (`AF_Genewiz_all_counts_no246.csv`, contrast Fetal-Adult, axes ±6/0–6) or cell counts (`AF_Cells_raw_counts.csv`, contrast Dual-Adult or Single-Adult, axes ±11/0–150). Currently configured for LCM.
- **Outputs:** `DEAnalysis_<group>_<contrast>.csv`, `SigUp_…` / `SigDown_…` CSVs (lines 257–265). The `ggsave` line is commented out, so the plot renders to the viewer only.
- **Dependencies / Notes:** Same as `Volcano_2.0.R`, plus a `plotBCV(z)` diagnostic. No existing GenericCode counterpart.

## `Volcano_Callouts_Significance.R`
- **Purpose:** "Overlay" volcano — highlight a custom gene list on the DE volcano (colored dots on a gray background) and export the filtered list.
- **Inputs:** LCM counts/metadata as above; callout list from `inputList = "EMT_Genes_Upregulated.csv"` in the input directory (line 29), `calloutList = "LIST"` (alternatives: `"TOP"` for top-50 DE genes, or `"GO"` for a GO-term gene set; lines 101–113).
- **Outputs** (`AF_Genewiz/Volcanoes/`, lines 330–340): `<inputContrast>_<listType>.png`, `<inputContrast>_Fetal-Adult...`-style `DEAnalysis_<contrast>.csv` of the full table, and `<nameOfOutputFile>_filtered_genes.csv` (callout genes passing FDR < 0.05 and |logFC| > `userLog2fc`).
- **Dependencies:** Same as `Volcano_2.0.R` plus `ggfortify`.
- **Usage/Notes:**
  - `nameOfPathway = "SYMBOL"` (line 38) chooses the key used to match the list; `"ENSEMBL"` also supported.
  - `colorStyle = "Custom"` colors all highlighted genes one color, auto-picked from the list name: red for `SigUp*`, blue for `SigDown*`, magenta otherwise (lines 57–64) — the current default input (`EMT_Genes_Upregulated.csv`) gives magenta. Switch to `"RedBlue"` for red-up/blue-down coloring.

## `Volcano_Callouts_Significance_Up&Down.R`
- **Purpose:** "Overlay" volcano comparing **cell** RNA-seq DE against **LCM-derived** gene lists — asks whether LCM-significant genes agree in direction in the cell data.
- **Inputs:** Cell counts `AF_Cells_raw_counts.csv` + `AF_Cells_metadata.csv` (lines 17–18), contrast `group: LAP-Adult` (lines 46–47); LCM lists `SigUp_DEAnalysis_adultfetal_Fetal-Adult.csv` / `SigDown_DEAnalysis_adultfetal_Fetal-Adult.csv` (lines 24–25; `listType = "SigUpandDown_LCM"`).
- **Outputs** (`AF_Genewiz/Volcanoes/`, lines 284–294): `LAP-Adult_SigUpandDown_LCM.png` and `LAP-Adult_SigUpandDown_LCM_filtered_genes.csv`. The plot title reports the overlap count vs total significant.
- **Dependencies:** Same as the other callout script.
- **Usage/Notes:**
  - Coloring is driven by the **colony** logFC direction (`ClogFC` joined in at the RedBlue branch): red = up in colonies, royalblue = down in colonies, gray = not in list, black = FDR > 0.05 — so a point's color tells you the LCM direction while its x-position shows the cell logFC.
  - The script's tail tallies concordance (`blueDown/blueUp/redDown/redUp` counts) for the console — the key QC numbers for this comparison.
  - Commented alternatives at the top switch to Dual-Adult / Single-Adult contrasts and their matching LCM lists.

## `Find_Exclusive_Genes.R`
- **Purpose:** Post-processing for the callout pipeline — intersect the Dual-Adult and Single-Adult filtered gene lists to find exclusive vs shared significant genes, and attach expression data for IPA.
- **Inputs:** `Dual-Adult_SigUp/SigDown_filtered_genes.csv` and `Single-Adult_SigUp/SigDown_filtered_genes.csv` in `AF_Genewiz/Volcanoes/` (lines 33–36), plus `DEAnalysis_group_Dual-Adult.csv` / `DEAnalysis_group_Single-Adult.csv` (lines 64–65).
- **Outputs** (Volcanoes dir for the set operations, input dir for the annotated ones; lines 52–95):
  - `Exclusive_Upregulated_Genes.csv`, `Exclusive_Downregulated_Genes.csv` (genes significant in only one comparison),
  - `Overlap_Upregulated_Genes.csv`, `Overlap_Downregulated_Genes.csv`,
  - `Dual_Exclusive_Upregulated/Downregulated_Genes.csv`, `Single_Exclusive_Upregulated/Downregulated_Genes.csv` (with logFC/PValue joined for IPA),
  - `All_Dual_Significant_Genes.csv`, `All_Single_Significant_Genes.csv`.
- **Dependencies:** `tidyverse` only.
- **Usage/Notes:** Run **after** producing the `*_filtered_genes.csv` files with the callout scripts (Dual-Adult and Single-Adult runs). Exclusivity = present in one list's merge but `NA` in the other (`allUpNA`/`allDownNA`, lines 44–45).

## Duplicate of a GenericCode version
- `Volcano from UAMS Proteomics (AF).R` — identical to `GenericCode/R Scripts/Volcano/Volcano from UAMS Proteomics (AF).R` — see that README.
