# Cluster figure scripts — CUT&Tag resubmission package

Ready-to-run on the cluster in the `cutntag` conda environment.
They regenerate every CUT&Tag figure panel from the pipeline outputs so
the figures match the (v2, final-run) manuscript text exactly, and they
produce the one panel that could not be made locally (Fig. 6B).

## Before running — confirm these paths

Both scripts take the pipeline output root from the `CNT_OUT`
environment variable or the default at the top of each file:

- `CNT_OUT` — default `~/BrownRIH/cnt_output`
  (the `OUTDIR` in the pipeline's `config.sh`). If you cleared and
  restored the outputs elsewhere, point `CNT_OUT` at the restored tree:
  `export CNT_OUT=/path/to/cnt_output`.
- Script 02 only: `GTF` — the Ensembl release 110 mRatBN7.2 GTF in the
  pipeline genome directory; the default filename is a guess at the
  exact name — check `ls $GENOME_DIR` and edit the one line if needed.
- Script 01 only: the Fig. 5A counts are read from
  `deseq2/<mark>/cells_Adult_vs_<X>.results.tsv`. If the per-population
  filenames differ, the script stops with a MISSING INPUT message
  naming the file — adjust `CONTRASTS_5A` / the filename pattern at the
  top.

## Run order

**One-shot (recommended):** copy `run_all_figures.sh` next to the two
Python scripts and run it. It restores the gene annotation file (GTF)
into the genome folder if missing (annotation only — no reference
FASTA or index; nothing is re-aligned), installs pyBigWig into the
cutntag environment, then runs script 01 and script 02 in order:

```bash
bash run_all_figures.sh
```

Folder locations default to the original cluster layout and can be
overridden: `CNT_OUT=/path/to/cnt_output GENOME_DIR=/path/to/genome bash run_all_figures.sh`.

Or run the steps manually:

```bash
conda activate cutntag
cd <this directory>
python 01_regen_panels.py          # no extra dependencies
pip install pyBigWig               # once, for script 02
python 02_fig6b_browser_tracks.py
```

Outputs land in `$CNT_OUT/figures_final/` (PDF + PNG each).

## What each script produces

**01_regen_panels.py** — all count/boxplot/scatter/summary panels:
- `fig5A_sig_peak_counts` — FDR<0.05 peak counts per adult-vs-fetal
  contrast, counted directly from the DESeq2 results tables.
- `fig5B_promoter_boxplots` — manuscript Fig. 5B (cells promoter
  log2FC by DEG class).
- `fig6A_retention_scatters` — manuscript Fig. 6A.
- `figA_promoter_medians` — concept Fig. A (cells vs tissue promoter
  medians). **P-values are annotated from the pipeline's own
  `geneset_tests.tsv`, never recomputed** — Python-recomputed Wilcoxon
  p-values disagree with R for the small-n marks (Oct 16 resolution).
- `figB_h3k36me1_signflip` — concept Fig. B (needs
  `sensitivity/3m1/tissue_nodedup.results.tsv`).
- `figC_qc_cells_vs_tissue` — concept Fig. C (from `qc/metrics.tsv`).
- `figD_concordance_summary` — concept Fig. D (deduplicated values
  computed from the retention tables; duplicate-retained values read
  from `sensitivity/concordance_summary.tsv`).

**02_fig6b_browser_tracks.py** — manuscript Fig. 6B:
`fig6b_tracks_cacna1b` (H3K27me1) and `fig6b_tracks_cacna1e`
(H3K27me3). Population-mean CPM tracks (Adult / LAP / Dual / Single,
mean of two replicates, 100 bp bins) over each gene ± 10 kb, with a
gene-model track. Gene coordinates are looked up in the GTF at
runtime. To change loci, populations, or flanks, edit the parameter
block at the top.

## Revised Oct 22, 2026

Both scripts were revised (headers list the changes): script 02 locks
one shared y-scale per locus for the browser tracks and uses the
manuscript mark names in titles; script 01 annotates the Fig. 5B
boxplot titles with the pipeline p-values (from geneset_tests.tsv,
as figA does), fixes the figA row collision with a taller, spaced
grid, and saves figB with a tight bounding box so its y-axis label
is no longer clipped. `run_all_figures.sh` is unchanged and still
drives both. Outputs overwrite the same filenames in
`$CNT_OUT/figures_final/`.

**Minimal restore list** (if `cnt_output` was cleared from the
cluster, these are the only pieces the figure scripts read):
- `cnt_output/integration/` (the four `*.retention_signature.tsv`
  + `geneset_tests.tsv`)
- `cnt_output/deseq2/` (per mark: the four
  `cells_Adult_vs_*.results.tsv` contrasts)
- `cnt_output/qc/metrics.tsv`
- `cnt_output/sensitivity/` (`concordance_summary.tsv` +
  `3m1/tissue_nodedup.results.tsv`)
- `cnt_output/bigwig/` — only the 16 isolated-cell files for marks
  2m1 and 2m3 (A1/A2/L1/L2/D1/D2/S1/S2 × 2m1/2m3) are read
- the genome GTF only (`Rattus_norvegicus.mRatBN7.2.110.gtf`) in
  the genome folder — no FASTA or index; `run_all_figures.sh`
  re-downloads it from Ensembl if it is missing.

## Notes

- These scripts supersede the local draft panels in `figure-drafts/`
  (which were plotted from uploaded tables; their Fig. A p-value
  annotations for two marks are known to be the Python-recomputed ones
  — do not use those annotations).
- Nothing here re-runs the pipeline; if `cnt_output` was cleared from
  the cluster, restore the copied folder first (or re-run the
  pipeline) — the scripts only read its outputs.
