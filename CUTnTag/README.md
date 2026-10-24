# CUT&Tag reanalysis — fetal hepatocyte transplant study

End-to-end reanalysis of the Azenta CUT&Tag data (48 samples: 32 pre-transplant
cell pellets + 16 post-transplant LCM tissue samples; 4 histone marks × 2 bioreps),
remapped to the current rat assembly (mRatBN7.2 / rn7).

## Why redo this

Azenta's delivered analysis (Jul–Aug 2023) was misconfigured: it merged peaks
*across different histone marks* per biological sample, ran replicate-vs-replicate
comparisons, and the key tissue comparison (adult host vs fetal colony, TA-vs-TF)
failed outright in their pipeline. The "0 differential peaks" result in the 2024
resubmission reflects that broken pipeline, not necessarily absent biology.

This pipeline replaces it with the standard modern CUT&Tag workflow:
consensus peak atlas **per mark** → fragment counting → DESeq2 differential
binding → promoter-level integration with the RNA-seq DEGs → the
**retention-signature test** (is the pre-transplant fetal-vs-adult mark landscape
preserved in post-transplant colony-vs-host tissue?).

## The biological hypothesis under test (from thesis Ch. 5)

PRC2/EZH2 deposits H3K27me3 (repressive, adult-high); H3K27me1 marks active
gene bodies (fetal-high); H3K36me2/3 antagonizes H3K27me3 (fetal-high).
Prediction: at promoters of polycomb-target / ion-transport genes
(CACNA1B, KCNQ3, CACNA1E, SLC6A1, SLC24A4, SLC6A5, ROBO3, GRIK1, …),
fetal colonies show **depleted H3K27me3** and **enriched H3K27me1/H3K36me2**
relative to adult host — i.e. the fetal mark landscape is retained 10 months
post-transplant.

## Layout

```
config.sh                 # <-- EDIT THIS: paths, threads, spike-in flag, DEG list
samples.tsv               # 48-sample manifest (fill in `pool` after checking SampleSheets)
envs/environment.yml      # conda env
scripts/00_setup.sh       # download rn7 + GTF, bowtie2 index, SEACR, chrom.sizes
scripts/01_qc.sh          # fastqc + fastp trim + multiqc
scripts/01b_poolmap.sh      # rebuild pool/lane map from FASTQ headers (Azenta never sent one)
scripts/02_align.sh       # bowtie2 → filter → markdup → dedup BAMs (+ optional spike-in)
scripts/03_peaks.sh       # SEACR peak calls, bigWigs, FRiP, QC metrics table
scripts/04_consensus_deseq2.R  # per-mark consensus peaks, counts, DESeq2 contrasts
scripts/05_integration.R  # annotate, DEG overlap, retention-signature test
scripts/06_metagenes.sh      # deeptools TSS profiles per mark × gene set
scripts/07_motifs.sh         # HOMER motif enrichment on fetal-retained peaks
run_all.sh                # driver
```

## Quick start

```bash
mamba env create -f envs/environment.yml && conda activate cutntag
# 1. edit config.sh (FASTQ_DIR, OUTDIR, GENOME_DIR, USE_SPIKEIN, DEG_CSV)
# 2. verify the `pool` column in samples.tsv against your SampleSheets
source config.sh
bash scripts/00_setup.sh
bash run_all.sh
```

## Key design decisions

- **Reference:** mRatBN7.2 (rn7), Ensembl release 110 — the original used Rnor6.0.
- **No IgG/input controls exist** (Azenta was told H3K27me3 is the positive
  control), so SEACR runs in no-control stringent mode; no input-subtraction.
- **Consensus peaks per mark, not replicate intersection.** Union of SEACR
  peaks within a mark → merged → kept if seen in ≥2 samples. Fragment counts
  per consensus peak → DESeq2. This is what rescues the tissue comparison.
- **Batch:** cell biorep 1 vs biorep 2 were sequenced as separate 16-sample
  pools (index conflicts), so biorep ≡ sequencing batch for cells. The design
  includes `pool` wherever a contrast spans pools. Tissue ran as its own
  project/run. Verify `pool` in samples.tsv against the SampleSheets.
- **Normalization:** no spike-in was used (confirmed Oct 2026 — Active Motif's
  spike-in normalization is incompatible with the CUT&Tag-IT kit), so DESeq2
  median-of-ratios is the primary normalization. Caveat: the mass-spec data
  show *global* mark-abundance shifts (e.g. H3K27me1 fetal-high), which violates
  that assumption — Azenta's own cell-vs-cell diffs show the artifact (18,118 vs
  266 peaks for Adult-vs-Fetal H3K27me1). Interpret differential-binding
  *directions* jointly with the MS abundances, not in isolation.
- **n=2 per group:** DESeq2 runs, but dispersion estimates are noisy. Report
  FDR<0.05 hits *and* the ranked list; lean on gene-set-level inference
  (retention signature, metagene plots) rather than individual peak calls.
- **Tissue QC gate:** samples with catastrophically low unique-fragment counts
  are flagged in `qc/metrics.tsv` and excluded from peak-level DESeq2, but
  retained (flagged) in aggregate plots.

## Contrasts run per mark (04)

Cells: Adult vs Fetal-pooled (L+D+S), Adult vs LAP, Adult vs Dual,
Adult vs Single, LAP vs Dual, LAP vs Single, Dual vs Single.
Tissue: TA (adult host) vs TF (fetal colony).
Retention: FetalCells vs TF, AdultCells vs TA (interpret cautiously —
different kits: cells vs tissue CUT&Tag-IT).

## Outputs

- `qc/multiqc_report.html`, `qc/metrics.tsv` (mapping/dup/FRiP/spike-in per sample)
- `peaks/*.stringent.bed`, `bigwig/*.bw` (CPM or spike-in normalized)
- `deseq2/<mark>/<contrast>.results.tsv`, PCA + volcano plots
- `integration/retention_signature.tsv`, metagene inputs, annotated peak tables
- `motifs/` HOMER results on fetal-retained promoter peaks

## What the Azenta reports tell us (parsed Oct 2026 — no recompute needed)

- Tissue (16 samples): 8–55M read pairs each, ~90% align, but **~90–97%
  duplication** → only ~0.6–1.5M usable pairs/sample. Per-sample peaks are
  real (300–2,700), but replicate overlap is ~5% (e.g. TF-2m3: 2,729+2,750 →
  123 merged) — this is exactly why consensus-then-test died and
  union-consensus + count-based DESeq2 is the fix. TA2-3m2 is the worst
  tissue sample (30.5% alignment).
- Cells: H3K27me1/me3 mostly fine (3–15M usable); H3K36me2 heavily
  duplicated — L2-3m2 (19K aligned) and S2-3m2 (157K usable) are failed
  samples, L1-3m2/S1-3m2/D2-3m2 marginal (~0.8–1.2M). Expect H3K36me2
  contrasts involving those samples to be weak — consistent with the
  manuscript's H3K36me2 heterogeneity note.
- Mitochondrial content negligible (~0.01%) — no mito filtering needed
  beyond the default.

## Still to fill in

1. Sanity-check `qc/pool_map.tsv` after 01b runs (expect ~3 runs:
   two 16-sample cell pools + one 16-sample tissue run).

DEG input is ready: `lcm_degs.csv` (953 rows: 786 colony-up + 167 host-up,
matching the thesis/manuscript counts exactly; symbols borrowed from the
annotated tables where available, joined to peaks by stable ENSEMBL ID).
Note: the annotated tables were built from a slightly different DE run
(860/186 rows, marginally different logFCs) — the 786/167 sets are the
manuscript-exact ones and are what's used.
