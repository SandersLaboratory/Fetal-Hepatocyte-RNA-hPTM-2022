# GO — README

Gene Ontology summarization and treemap figures for the Adult–Fetal project. Two approaches are used here: the `rrvgo` R package (re-runnable from GO output tables) and frozen REVIGO-server exports (fixed data baked into the script).

## `RRVGO.R`
- **Purpose:** Summarize and de-redundantize long GO-enrichment lists using the `rrvgo` package (semantic-similarity clustering), then plot the reduced set as a treemap (plus optional heatmap/scatter/wordcloud).
- **Inputs:** `Colony GO Output - pval 10^-3.txt` and `Host GO Output - pval 10^-3.txt` (tab-delimited, with `GO.term` and `pValue` columns) in `Adult-Fetal/Data Files/go/`.
- **Outputs:** Plots print to the RStudio viewer (no files written).
- **Dependencies:** `tidyverse`, `rrvgo`, `org.Rn.eg.db` (used via the `organism` setting; Bioconductor).
- **Usage/Notes:** Configure the "User Settings" block (lines 8–27):
  - `whichSet = "host"` (or `"colony"`) selects which file is summarized.
  - `pvalCutoff = 1e-5` filters terms; `GOthreshold = 0.7` sets the semantic-similarity reduction threshold; `ont = "BP"` (hard-coded at line ~50) restricts to Biological Process; `method = "Rel"`.
  - `whichPlot = "treemap"` (or `"all"`, `"scatter"`, `"heatmap"`, `"wordcloud"`).
  - `organism` is currently `"org.Rn.eg.db"` (rat) — switch to `"org.Hs.eg.db"` for human gene sets.
  - Authored by Anders Ohman, Spring 2022 (header comment).

## REVIGO treemap exports
These four scripts are exports from the REVIGO web server (revigo.irb.hr; credit header inside each file, citing Supek et al. PLoS ONE 2011), with the GO term table hard-coded as `revigo.data`. They render one static treemap PDF in the working directory via the `treemap` package. Because the data is baked in, they are **not re-runnable on new GO results** — for new data, re-export a fresh script from REVIGO (or use `RRVGO.R` above).

- **`Host_RevigoTreeMap.R`** → `host_revigo_treemap.pdf` (line 46). 14 host GO terms, dominated by metabolic/catabolic processes. Note: this directory also contains rendered outputs (`host_revigo_treemap.pdf/.png`, `host_revigo_treemap1/2.pdf/.png`) from earlier runs.
- **`Colony_RevigoTreeMap.R`** → `colony_revigo_treemap.pdf` (line 53). 21 colony GO terms, dominated by ion/transmembrane transport. Also rendered in this directory (`colony_revigo_treemap.pdf/.png`, `colony_revigo_treemap1.pdf/.png`).
- **`Metabolism - RevigoTreeMap.R`** → `metabolism_revigo_treemap.pdf` (line 48). 16 metabolism-related GO terms (host-flavored: carboxylic/lipid acid metabolism). Note the plot **title still says "Host GO Processes"** — leftover from the copied export.
- **`Ion Transporter - RevigoTreeMap.R`** → `iontransport_revigo_treemap.pdf` (line 52). 20 ion-transport GO terms (colony-flavored). Note the plot **title still says "Colony GO Processes"** — leftover from the copied export.
