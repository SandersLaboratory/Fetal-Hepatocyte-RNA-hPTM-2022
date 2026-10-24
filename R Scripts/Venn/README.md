# Venn — README

Venn diagrams (via `VennDiagram`) comparing significant gene lists across AF datasets.

## `Generic Venn Diagram Function.R`
- **Purpose:** Copy-paste template for 2-, 3-, or 4-way Venn diagrams. Fill in `number_of_categories`, the set labels (`set1`…`set4`), and the item vectors (`s1`…`s4`); the script draws the matching diagram (lines 19–32).
- **Inputs:** Hard-coded example vectors at the top (lines 11–16) — replace with your gene lists.
- **Outputs:** Plot renders to the graphics device (no file saved).
- **Dependencies:** `VennDiagram`.
- **Usage/Notes:** Counts and all pairwise/higher-order intersections are computed from the vectors with `length()`/`intersect()`, so areas are exact for the input lists. Differs from the GenericCode version (same filename, different size — the example vectors differ), so it is documented here rather than by cross-reference.

## `Venn from Microarray vs RNA-Seq.R`
- **Purpose:** Compare significant (FDR < 0.05) DE gene lists from microarray (Dual = LAP+/OC2+, Single = LAP+/OC2−) against RNA-seq (Colony): a 3-way Venn plus a Dual-vs-Colony 2-way Venn, and export the overlap gene sets.
- **Inputs:**
  - `Adult_Dual/Adult_Dual_annotated_DE_Fetal-Adult.csv` and `Adult_Single/Adult_Single_annotated_DE_Fetal-Adult.csv` in `Adult-Fetal/Data Files/Microarray/` (lines ~19–22; filtered on `adj.P.Val`).
  - `DEAnalysis_adultfetal_Fetal-Adult.csv` in `Adult-Fetal/Data Files/rna-seq/AF_Genewiz/` (line 24; ENSEMBL keys annotated to SYMBOL + UNIPROT via `org.Rn.eg.db`, filtered on `FDR`, at lines 26–38).
- **Outputs** (`Adult-Fetal/Data Files/venn/`, lines 79–85): `Genes_In_Dual&Single&Colony.csv`, `Genes_In_Dual&Colony.csv`, `Genes_In_Single&Colony.csv`, `Genes_Only_In_Dual&Colony_NotSingle.csv`.
- **Dependencies:** `VennDiagram`, `tidyverse`, `heatmaply`, `org.Hs.eg.db`, `org.Rn.eg.db`, `annotate`, `rstudioapi`.
- **Usage/Notes:** NA symbols are dropped before comparison (lines 50–53); colony symbols come from the `genesymbol` column, so exact column naming matters. Differs from the GenericCode version (same filename, different size).

## Duplicates of GenericCode versions
- `Venn from serum Proteomics.R` — identical to `GenericCode/R Scripts/Venn/Venn from serum Proteomics.R` — see that README.
- `Venn from Proteomics.R` — identical to `GenericCode/R Scripts/Venn/Venn from Proteomics.R` — see that README.
