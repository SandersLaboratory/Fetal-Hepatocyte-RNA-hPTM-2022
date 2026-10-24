#!/usr/bin/env Rscript
# 08_sensitivity_tissue_nodedup.R — companion to 08_sensitivity_tissue_nodedup.sh.
# Recounts each mark's EXISTING consensus atlas on the duplicate-retaining
# tissue BAMs, re-runs TA vs TF in DESeq2, and recomputes concordance with
# the (unchanged) cells Adult-vs-FetalPooled effects. Writes per-mark
# results + one summary table under $OUTDIR/sensitivity/.
suppressPackageStartupMessages({
  library(Rsubread); library(DESeq2); library(readr); library(dplyr)
})
SAMPLES_TSV <- Sys.getenv("SAMPLES_TSV"); OUTDIR <- Sys.getenv("OUTDIR")
THREADS <- as.integer(Sys.getenv("THREADS"))
samp <- read_tsv(SAMPLES_TSV, show_col_types = FALSE)
SDIR <- file.path(OUTDIR, "sensitivity"); dir.create(SDIR, recursive = TRUE, showWarnings = FALSE)

summ <- list()
for (mk in c("2m1", "2m3", "3m1", "3m2")) {
  message("==== mark ", mk, " ====")
  bed <- read_tsv(file.path(OUTDIR, "deseq2", mk, "consensus_peaks.bed"),
                  col_names = c("chr", "start", "end"), show_col_types = FALSE)
  saf <- data.frame(GeneID = paste0("peak", seq_len(nrow(bed))),
                    Chr = bed$chr, Start = bed$start + 1L, End = bed$end, Strand = "*")
  ss <- samp %>% filter(mark == mk, type == "tissue")
  bams <- file.path(SDIR, "bam", paste0(ss$sample, ".nodedup.bam"))
  fc <- featureCounts(bams, annot.ext = saf, isGTFAnnotationFile = FALSE,
                      isPairedEnd = TRUE, countReadPairs = TRUE, nthreads = THREADS)
  counts <- fc$counts; colnames(counts) <- ss$sample
  message("total fragments counted (nodedup): ", sum(counts))

  ss$cond <- factor(ss$population, levels = c("TissueAdult", "TissueFetal"))
  cd <- DataFrame(ss); rownames(cd) <- ss$sample
  dds <- DESeq(DESeqDataSetFromMatrix(counts, cd, ~ cond), quiet = TRUE)
  res <- results(dds, contrast = c("cond", "TissueFetal", "TissueAdult"))  # + = colony-high
  out <- data.frame(chr = bed$chr, start = bed$start, end = bed$end,
                    baseMean = res$baseMean, log2FC_mle = res$log2FoldChange,
                    pvalue = res$pvalue, padj = res$padj)
  mdir <- file.path(SDIR, mk); dir.create(mdir, showWarnings = FALSE)
  write_tsv(out, file.path(mdir, "tissue_nodedup.results.tsv"))
  message("nodedup TA vs TF: ", sum(out$padj < 0.05, na.rm = TRUE), " FDR<0.05 peaks")

  # concordance with the cells contrast (same consensus atlas; join on coords)
  rc <- read_tsv(file.path(OUTDIR, "deseq2", mk, "cells_Adult_vs_FetalPooled.results.tsv"),
                 show_col_types = FALSE)
  cells_lfc <- if ("log2FC_mle" %in% names(rc)) rc$log2FC_mle else rc$log2FC
  j <- data.frame(cells = cells_lfc, tissue = out$log2FC_mle) %>%
    filter(is.finite(cells), is.finite(tissue))
  agree <- mean(sign(j$cells) == sign(j$tissue))
  strong <- j %>% filter(abs(cells) > 1)
  agree_s <- mean(sign(strong$cells) == sign(strong$tissue))
  rho <- suppressWarnings(cor(j$cells, j$tissue, method = "spearman"))

  # DEG-promoter gene-set test on the nodedup tissue effects, reusing the
  # annotation/deg_class from step 05's retention table (same atlas order
  # is NOT assumed — join on coordinates)
  rt <- read_tsv(file.path(OUTDIR, "integration", paste0(mk, ".retention_signature.tsv")),
                 show_col_types = FALSE)
  m <- rt %>% select(chr, start, end, annotation, deg_class) %>%
    inner_join(out %>% select(chr, start, end, tissue_nd = log2FC_mle),
               by = c("chr", "start", "end")) %>%
    filter(grepl("^Promoter", annotation), !is.na(deg_class))
  wp <- NA_real_
  if (sum(m$deg_class == "colony_up") > 2 && sum(m$deg_class != "colony_up") > 2)
    wp <- wilcox.test(m$tissue_nd[m$deg_class == "colony_up"],
                      m$tissue_nd[m$deg_class != "colony_up"])$p.value
  summ[[mk]] <- data.frame(mark = mk, n_peaks = nrow(j), sign_agreement = agree,
                           sign_agreement_strong_cells = agree_s, spearman_rho = rho,
                           promoter_colonyup_wilcox_p = wp)
  message(sprintf("concordance: sign agreement %.3f (strong-cells %.3f), Spearman %.3f, promoter colony_up Wilcox p %s",
                  agree, agree_s, rho, ifelse(is.na(wp), "NA", sprintf("%.3f", wp))))
}
write_tsv(bind_rows(summ), file.path(SDIR, "concordance_summary.tsv"))
message("sensitivity summary written to ", file.path(SDIR, "concordance_summary.tsv"))
