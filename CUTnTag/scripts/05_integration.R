#!/usr/bin/env Rscript
# 05_integration.R — annotate peaks, join RNA-seq DEGs, retention-signature test
#
# Peak annotation is done directly with GenomicRanges (no ChIPseeker: its
# dependency subtree pulls in system font libraries and fragile graphics
# packages that fail to build on this cluster):
#   * TSSs parsed from Ensembl GTF "gene" features (1-bp, strand-aware)
#   * each consensus peak assigned to its nearest TSS (either strand)
#   * signed distanceToTSS = peak midpoint minus TSS position, sign flipped on
#     minus-strand genes (negative = upstream, positive = downstream)
#   * promoter = |distance| <= 3 kb, binned <=1kb / 1-2kb / 2-3kb
#   * non-promoter peaks labelled Intragenic (overlaps a gene body) or Intergenic
suppressPackageStartupMessages({
  library(GenomicRanges); library(readr); library(dplyr); library(ggplot2)
})
OUTDIR <- Sys.getenv("OUTDIR"); GENOME_DIR <- Sys.getenv("GENOME_DIR")
DEG_CSV <- Sys.getenv("DEG_CSV")
MARKS <- c("2m1" = "H3K27me1", "2m3" = "H3K27me3", "3m1" = "H3K36me1", "3m2" = "H3K36me2")
idir <- file.path(OUTDIR, "integration"); dir.create(idir, showWarnings = FALSE)
gsdir <- file.path(idir, "genesets"); dir.create(gsdir, showWarnings = FALSE)
hdir <- file.path(idir, "homer_beds"); dir.create(hdir, showWarnings = FALSE)

# --- gene models straight from the GTF used for alignment ---
gtf <- list.files(GENOME_DIR, pattern = "mRatBN7.2.*\\.gtf$", full.names = TRUE)[1]
gff <- read_tsv(gtf, comment = "#", col_names = FALSE, show_col_types = FALSE,
                col_types = "ccciiiccc")  # seqname source feature start end score strand frame attr
g <- gff[gff$X3 == "gene", ]
g$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", g$X9)
g$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", g$X9)
g$gene_name[g$gene_name == g$X9] <- NA  # no gene_name attribute -> NA
g <- g[!duplicated(g$gene_id), ]
genes <- data.frame(gene_id = g$gene_id, gene_name = g$gene_name)

tss_pos <- ifelse(g$X7 == "+", g$X4, g$X5)  # GTF is 1-based closed
tss <- GRanges(seqnames = g$X1, ranges = IRanges(tss_pos, tss_pos), strand = g$X7,
               gene_id = g$gene_id, gene_name = g$gene_name)
genebody <- GRanges(seqnames = g$X1, ranges = IRanges(g$X4, g$X5), strand = g$X7)
message("genes in GTF: ", length(tss))

annotate_peaks <- function(peaks) {
  # distanceToNearest returns NO hit for peaks on scaffolds without any
  # annotated gene (unplaced JACYVU/MU contigs in the toplevel assembly:
  # 140 of 16,730 H3K27me1 peaks). Pre-fill one NA row per peak and fill in
  # by query index, so the result always aligns row-for-row with `peaks`.
  out <- data.frame(geneId = rep(NA_character_, length(peaks)),
                    symbol = rep(NA_character_, length(peaks)),
                    annotation = rep(NA_character_, length(peaks)),
                    distanceToTSS = rep(NA_integer_, length(peaks)))
  nn <- distanceToNearest(peaks, tss)
  qi <- queryHits(nn); si <- subjectHits(nn)
  th <- tss[si]
  mid <- (start(peaks)[qi] + end(peaks)[qi]) %/% 2L
  tpos <- start(th)
  signed <- ifelse(as.character(strand(th)) == "+", mid - tpos, tpos - mid)
  ad <- abs(signed)
  annotation <- ifelse(ad <= 1000, "Promoter (<=1kb)",
                ifelse(ad <= 2000, "Promoter (1-2kb)",
                ifelse(ad <= 3000, "Promoter (2-3kb)", "Distal")))
  intra <- overlapsAny(peaks[qi], genebody)
  annotation[annotation == "Distal"] <-
    ifelse(intra[annotation == "Distal"], "Intragenic", "Intergenic")
  out$geneId[qi] <- th$gene_id
  out$symbol[qi] <- th$gene_name
  out$annotation[qi] <- annotation
  out$distanceToTSS[qi] <- signed
  out
}

# --- RNA-seq DEGs (LCM colony vs host) ---
deg <- read_csv(DEG_CSV, show_col_types = FALSE) %>%
  mutate(deg_class = ifelse(FDR < 0.05 & logFC > 0, "colony_up",
                     ifelse(FDR < 0.05 & logFC < 0, "host_up", "ns")))

LEADING_EDGE <- c("CACNA1B", "KCNQ3", "CACNA1E", "SLC6A1", "SLC24A4", "SLC6A5",
                  "ROBO3", "GRIK1")  # thesis Ch.5 polycomb GSEA leading edge
NOMINATED <- c("CACNA1B", "KCNJ12", "SLC34A3", "TMC7", "WNT3A")  # thesis validation targets

# --- gene-set TSS BEDs for deeptools metagenes ---
write_tss_bed <- function(symbols, file) {
  gid <- genes$gene_id[match(toupper(symbols), toupper(genes$gene_name))]
  tt <- tss[tss$gene_id %in% na.omit(gid)]
  bed <- data.frame(chr = as.character(seqnames(tt)), start = start(tt) - 1L,
                    end = end(tt), name = tt$gene_id, score = 0,
                    strand = as.character(strand(tt)))
  write_tsv(bed, file, col_names = FALSE)
}
write_tss_bed(deg$gene_symbol[deg$deg_class == "colony_up"], file.path(gsdir, "colony_up.tss.bed"))
write_tss_bed(deg$gene_symbol[deg$deg_class == "host_up"],   file.path(gsdir, "host_up.tss.bed"))
write_tss_bed(LEADING_EDGE, file.path(gsdir, "leading_edge.tss.bed"))
write_tss_bed(genes$gene_name, file.path(gsdir, "background.tss.bed"))

tests <- list(); retention_all <- list()
for (mk in names(MARKS)) {
  message("==== ", MARKS[[mk]], " ====")
  mdir <- file.path(OUTDIR, "deseq2", mk)
  bed <- read_tsv(file.path(mdir, "consensus_peaks.bed"), col_names = FALSE, show_col_types = FALSE)
  peaks <- GRanges(bed$X1, IRanges(bed$X2 + 1L, bed$X3))
  ann <- annotate_peaks(peaks)

  rc <- read_tsv(file.path(mdir, "cells_Adult_vs_FetalPooled.results.tsv"), show_col_types = FALSE)  # + = fetal-high
  rt <- read_tsv(file.path(mdir, "tissue_TissueAdult_vs_TissueFetal.results.tsv"), show_col_types = FALSE)  # + = colony-high
  # Direction/concordance analyses use UNSHRUNKEN (MLE) log2FCs: apeglm
  # shrinkage is information-dependent and reduces the low-count tissue
  # contrast to ~0, which would erase its directions. Fall back to the
  # shrunken column only for result files written before log2FC_mle existed.
  lfc_of <- function(r) if ("log2FC_mle" %in% names(r)) r$log2FC_mle else r$log2FC
  tab <- data.frame(chr = rc$chr, start = rc$start, end = rc$end,
                    gene_id = ann$geneId, symbol = ann$symbol,
                    annotation = ann$annotation, distToTSS = ann$distanceToTSS,
                    cells_log2FC = lfc_of(rc), cells_padj = rc$padj,
                    tissue_log2FC = lfc_of(rt), tissue_padj = rt$padj) %>%
    left_join(deg %>% select(ensembl_id, deg_logFC = logFC, deg_FDR = FDR, deg_class),
              by = c("gene_id" = "ensembl_id"))
  tab$mark <- MARKS[[mk]]
  write_tsv(tab, file.path(idir, paste0(mk, ".retention_signature.tsv")))
  retention_all[[mk]] <- tab

  prom <- tab %>% filter(grepl("^Promoter", annotation))
  # scatter: pre-transplant fetal-vs-adult vs post-transplant colony-vs-host (promoter peaks)
  q <- ggplot(prom, aes(cells_log2FC, tissue_log2FC, color = deg_class)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = c(colony_up = "firebrick", host_up = "steelblue", ns = "grey70"),
                       na.value = "grey70") +
    labs(title = paste0(MARKS[[mk]], ": retention of fetal mark landscape at promoters"),
         x = "pre-transplant: fetal vs adult (log2FC)", y = "post-transplant: colony vs host (log2FC)") +
    theme_minimal()
  ggsave(file.path(idir, paste0(mk, ".retention_scatter.png")), q, width = 7, height = 6)

  # gene-set test: is the colony_up DEG set shifted in the fetal direction?
  for (cls in c("colony_up", "host_up")) {
    gg <- prom %>% filter(!is.na(deg_class))
    # explicit two-vector form: the estimate then reads in-set minus
    # out-of-set. (The formula form `y ~ (deg_class == cls)` silently
    # reports the reverse direction — first factor level minus second.)
    wt <- wilcox.test(gg$tissue_log2FC[gg$deg_class == cls],
                      gg$tissue_log2FC[gg$deg_class != cls], conf.int = TRUE)
    wc <- wilcox.test(gg$cells_log2FC[gg$deg_class == cls],
                      gg$cells_log2FC[gg$deg_class != cls], conf.int = TRUE)
    tests[[length(tests) + 1]] <- data.frame(
      mark = MARKS[[mk]], gene_set = cls, n_set = sum(gg$deg_class == cls),
      tissue_p = wt$p.value, tissue_median_shift = wt$estimate,
      cells_p = wc$p.value, cells_median_shift = wc$estimate)
  }

  # fetal-retained promoter peaks -> HOMER BEDs
  # permissive marks (fetal-high): concordant positive; repressive H3K27me3: concordant negative
  sgn <- if (mk == "2m3") -1 else 1
  fb <- prom %>% filter(deg_class == "colony_up",
                        sign(cells_log2FC) == sgn, sign(tissue_log2FC) == sgn,
                        pmin(cells_padj, tissue_padj, na.rm = TRUE) < 0.1)
  message("fetal-retained promoter peaks for HOMER: ", nrow(fb))
  write_tsv(data.frame(chr = fb$chr, start = fb$start, end = fb$end),
            file.path(hdir, paste0(mk, ".fetal_retained.bed")), col_names = FALSE)
}
write_tsv(bind_rows(tests), file.path(idir, "geneset_tests.tsv"))
write_tsv(bind_rows(retention_all) %>% filter(toupper(symbol) %in% c(LEADING_EDGE, NOMINATED)),
        file.path(idir, "leading_edge_nominated.tsv"))
message("05 done")
