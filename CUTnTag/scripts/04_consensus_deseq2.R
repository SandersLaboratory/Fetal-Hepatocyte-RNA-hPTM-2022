#!/usr/bin/env Rscript
# 04_consensus_deseq2.R — per-mark consensus peaks, fragment counts, DESeq2 contrasts
suppressPackageStartupMessages({
  library(GenomicRanges); library(Rsubread); library(DESeq2)
  library(ggplot2); library(pheatmap); library(readr); library(dplyr)
})
SAMPLES_TSV <- Sys.getenv("SAMPLES_TSV"); OUTDIR <- Sys.getenv("OUTDIR")
THREADS <- as.integer(Sys.getenv("THREADS")); USE_SPIKEIN <- Sys.getenv("USE_SPIKEIN") == "true"
MIN_FRAGS <- 100000   # QC gate: samples below this are excluded from DESeq2 (still in aggregate plots).
                      # 100K floor (Anders, Oct 13, 2026) drops only the dead L2-3m2 (42K frags;
                      # its FASTQ held just 35K reads) and keeps marginal-but-real S2-3m2 (350K).
GAP <- 200           # merge gap for broad histone peaks

samp <- read_tsv(SAMPLES_TSV, show_col_types = FALSE)
met  <- read_tsv(file.path(OUTDIR, "qc/metrics.tsv"), show_col_types = FALSE)
keep <- met$sample[met$unique_frags >= MIN_FRAGS]
drop <- setdiff(met$sample, keep)
if (length(drop)) message("QC gate dropped from DESeq2: ", paste(drop, collapse = ", "))
samp <- samp %>% filter(sample %in% keep)
bam_of <- function(s) file.path(OUTDIR, "bam", paste0(s, ".dedup.bam"))

# pool mapping: FASTQ headers (01b_poolmap.sh) override the samples.tsv guess,
# since Azenta never documented which 16+16 split they sequenced
pm <- file.path(OUTDIR, "qc/pool_map.tsv")
if (file.exists(pm)) {
  pmap <- read_tsv(pm, show_col_types = FALSE)
  samp$pool <- paste0("fc", pmap$flowcell, "_lane", pmap$lane)[match(samp$sample, pmap$sample)]
  message("pool assignments taken from FASTQ headers")
}

read_bed3 <- function(f) {
  # weak samples may have no peaks file (or an empty one) — return an empty
  # GRanges rather than aborting the run
  if (!file.exists(f)) return(GRanges())
  d <- read_tsv(f, col_names = FALSE, show_col_types = FALSE)
  if (nrow(d) == 0) return(GRanges())
  GRanges(d$X1, IRanges(d$X2 + 1L, d$X3))  # BED 0-based -> 1-based
}

# condition codings used across contrasts
samp <- samp %>% mutate(
  fetal_pooled = ifelse(population %in% c("LAP", "Dual", "Single"), "FetalPooled", population)
)
# No cross-project (cells-vs-tissue) DESeq2 contrasts: pool/project/sample-type
# are fully confounded (all tissue = one pool; each pool maps to one condition),
# so ~pool + cond is rank-deficient and DESeq2 refuses to fit it — and the
# comparison would be uninterpretable anyway. The retention question is tested
# in 05 by comparing within-project effect directions (cells vs tissue).
contrasts <- list(
  cells  = list(c("Adult", "FetalPooled"), c("Adult", "LAP"), c("Adult", "Dual"),
                c("Adult", "Single"), c("LAP", "Dual"), c("LAP", "Single"), c("Dual", "Single")),
  tissue = list(c("TissueAdult", "TissueFetal"))
)
# For "cells", per-population contrasts use `population`; only the pooled
# contrast (trt == "FetalPooled") uses the pooled coding. (Returning
# fetal_pooled for every cell contrast made the per-population filters match
# only the 2 Adult samples, so all of them were silently skipped as "too few".)
cond_col <- function(s, which, trt = NULL) {
  if (which == "cells" && !is.null(trt) && trt == "FetalPooled") return(s$fetal_pooled)
  s$population
}

for (mk in c("2m1", "2m3", "3m1", "3m2")) {
  message("==== mark ", mk, " ====")
  mdir <- file.path(OUTDIR, "deseq2", mk); dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
  ss <- samp %>% filter(mark == mk)
  grl <- lapply(ss$sample, function(s) read_bed3(file.path(OUTDIR, "peaks", paste0(s, ".seacr.stringent.bed"))))
  names(grl) <- ss$sample
  union_p <- reduce(unlist(GRangesList(grl)), min.gapwidth = GAP)
  if (length(union_p) == 0) { message("no SEACR peaks for mark ", mk, " — skipping mark"); next }
  present <- rowSums(sapply(grl, function(g) countOverlaps(union_p, g) > 0))
  consensus <- union_p[present >= 2]
  message("consensus peaks (>=2 samples): ", length(consensus))
  bed <- data.frame(chr = as.character(seqnames(consensus)),
                    start = start(consensus) - 1L, end = end(consensus))
  write_tsv(bed, file.path(mdir, "consensus_peaks.bed"), col_names = FALSE)

  saf <- data.frame(GeneID = paste0("peak", seq_along(consensus)),
                    Chr = bed$chr, Start = bed$start + 1L, End = bed$end, Strand = "*")
  fc <- featureCounts(bam_of(ss$sample), annot.ext = saf, isGTFAnnotationFile = FALSE,
                      isPairedEnd = TRUE, nthreads = THREADS)
  counts <- fc$counts; colnames(counts) <- ss$sample
  saveRDS(list(counts = counts, consensus = consensus, samples = ss),
          file.path(mdir, "counts.rds"))

  # PCA on vst (all samples of this mark)
  cd0 <- DataFrame(ss); rownames(cd0) <- ss$sample
  d0 <- DESeqDataSetFromMatrix(counts, cd0, ~1)
  vsd <- vst(d0, blind = TRUE)
  p <- plotPCA(vsd, intgroup = c("population", "type")) + ggtitle(paste0(mk, ": PCA"))
  ggsave(file.path(mdir, "pca.png"), p, width = 7, height = 5)

  for (grp in names(contrasts)) {
    for (cc in contrasts[[grp]]) {
      ctl <- cc[1]; trt <- cc[2]
      sub <- ss %>% filter(cond_col(., grp, trt) %in% c(ctl, trt))
      if (nrow(sub) < 3) { message("skip ", grp, " ", ctl, " vs ", trt, " (too few)"); next }
      sub$cond <- factor(cond_col(sub, grp, trt), levels = c(ctl, trt))  # log2FC = trt vs ctl
      # include pool only when the design stays estimable: ~pool + cond has
      # npool + 1 parameters, so require more samples than that (the 2-vs-1
      # H3K36me2 cell contrasts after the QC gate would otherwise saturate)
      des <- if (n_distinct(sub$pool) > 1 && nrow(sub) > n_distinct(sub$pool) + 1) ~pool + cond else ~cond
      cd <- DataFrame(sub); rownames(cd) <- sub$sample
      dds <- DESeqDataSetFromMatrix(counts[, sub$sample], cd, des)
      if (USE_SPIKEIN) {
        sp <- met$spikein_frac[match(sub$sample, met$sample)]
        sizeFactors(dds) <- sp / exp(mean(log(sp)))
      }
      dds <- DESeq(dds, quiet = TRUE)
      rn <- resultsNames(dds); coef <- rn[length(rn)]
      res <- as.data.frame(results(dds, name = coef))
      shr <- tryCatch(as.data.frame(lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)),
                      error = function(e) res)
      out <- data.frame(chr = bed$chr, start = bed$start, end = bed$end,
                        baseMean = res$baseMean,
                        # log2FC = apeglm-shrunken (reporting/volcanoes);
                        # log2FC_mle = unshrunken. Step 05's direction and
                        # concordance tests MUST use the MLE: shrinkage is
                        # information-dependent and shrinks the low-count
                        # tissue contrast to ~0 (observed max |LFC| ~1e-6),
                        # erasing its effect directions entirely.
                        log2FC = shr$log2FoldChange, log2FC_mle = res$log2FoldChange,
                        lfcSE = shr$lfcSE,
                        pvalue = res$pvalue, padj = res$padj)
      tag <- paste0(grp, "_", ctl, "_vs_", trt)
      write_tsv(out, file.path(mdir, paste0(tag, ".results.tsv")))
      # volcano
      out$sig <- ifelse(!is.na(out$padj) & out$padj < 0.05, "FDR<0.05",
                 ifelse(!is.na(out$pvalue) & out$pvalue < 0.01, "p<0.01", "ns"))
      v <- ggplot(out, aes(log2FC, -log10(pvalue), color = sig)) +
        geom_point(size = 0.6, alpha = 0.7) +
        scale_color_manual(values = c("FDR<0.05" = "red", "p<0.01" = "orange", "ns" = "grey60")) +
        ggtitle(paste0(mk, ": ", trt, " vs ", ctl)) + theme_minimal()
      ggsave(file.path(mdir, paste0(tag, ".volcano.png")), v, width = 6, height = 5)
      message(tag, ": ", sum(out$padj < 0.05, na.rm = TRUE), " FDR<0.05 peaks")
    }
  }
}
message("04 done")
