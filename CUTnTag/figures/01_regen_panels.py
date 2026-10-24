#!/usr/bin/env python3
"""
01_regen_panels.py — regenerate ALL CUT&Tag figure panels from cluster outputs.

Produces (PDF + PNG each), in OUT_FIGDIR:
  fig5A_sig_peak_counts      significant-peak counts per adult-vs-fetal contrast
  fig5B_promoter_boxplots    DEG-promoter log2FC distributions, cells (manuscript Fig. 5B)
  fig6A_retention_scatters   cells vs tissue log2FC scatter per mark (manuscript Fig. 6A)
  figA_promoter_medians      concept Fig. A: promoter medians, cells vs tissue, pipeline p-values
  figB_h3k36me1_signflip     concept Fig. B: H3K36me1 dedup vs duplicate-retained tissue medians
  figC_qc_cells_vs_tissue    concept Fig. C: per-sample FRiP / duplication / unique fragments
  figD_concordance_summary   concept Fig. D: sign agreement + Spearman, dedup vs dup-retained

IMPORTANT (per the Oct 16 v1/v2 resolution): p-values annotated on panels are
read from the pipeline's own geneset_tests.tsv (R wilcox.test) — never
recomputed here. Python-recomputed Wilcoxon p-values disagree with R for the
small-n marks (H3K27me1, H3K36me2); medians ARE computed from the retention
tables (they reproduce the pipeline values exactly).

REVISED Oct 22, 2026 — three fixes over the original version:
  1. fig5B boxplot titles now carry the pipeline p-value for each mark
     (same geneset_tests.tsv source as figA's annotations).
  2. figA (promoter medians): taller grid + explicit row spacing +
     tight_layout — the two rows' titles and tick labels collided in
     the original render.
  3. figB (sign-flip) is saved with bbox_inches="tight" — its y-axis
     label was clipped at the figure edge in the original render.
All other panels are unchanged.

Run inside the cutntag conda env on the cluster:
  conda activate cutntag
  python 01_regen_panels.py
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ PARAMETERS — CONFIRM BEFORE RUNNING ============
# Root of the pipeline output tree on the cluster (config.sh OUTDIR).
CNT_OUT   = os.environ.get("CNT_OUT", "~/BrownRIH/cnt_output")  # <-- CONFIRM
OUT_FIGDIR = os.path.join(CNT_OUT, "figures_final")                              # <-- CONFIRM or change
# Contrast result files under deseq2/<mark>/ used for the Fig. 5A counts.
# Verify the per-population filenames with:  ls $CNT_OUT/deseq2/2m1/
CONTRASTS_5A = ["FetalPooled", "LAP", "Dual", "Single"]   # file: cells_Adult_vs_<X>.results.tsv
SCATTER_MAX_POINTS = 8000        # per-mark subsample cap for fig6A (file size)
MARKS = [("2m1", "H3K27me1"), ("2m3", "H3K27me3"), ("3m1", "H3K36me1"), ("3m2", "H3K36me2")]
# =============================================================================

os.makedirs(OUT_FIGDIR, exist_ok=True)
INTEG = os.path.join(CNT_OUT, "integration")
DESQ  = os.path.join(CNT_OUT, "deseq2")
SENS  = os.path.join(CNT_OUT, "sensitivity")

def need(path):
    if not os.path.exists(path):
        sys.exit(f"MISSING INPUT: {path}\nCheck CNT_OUT and the pipeline output layout.")
    return path

ret = {}
for mk, _ in MARKS:
    ret[mk] = pd.read_csv(need(os.path.join(INTEG, f"{mk}.retention_signature.tsv")), sep="\t")
gst = pd.read_csv(need(os.path.join(INTEG, "geneset_tests.tsv")), sep="\t")
metrics = pd.read_csv(need(os.path.join(CNT_OUT, "qc", "metrics.tsv")), sep="\t")

def save(fig, name, tight=False):
    kw = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(os.path.join(OUT_FIGDIR, name + ".pdf"), **kw)
    fig.savefig(os.path.join(OUT_FIGDIR, name + ".png"), dpi=200, **kw)
    plt.close(fig)
    print("wrote", name)

def pipe_p(mark_label, gs, which):
    """Pipeline p-value from geneset_tests.tsv (never recomputed here)."""
    row = gst[(gst.mark == mark_label) & (gst.gene_set == gs)]
    assert len(row) == 1, f"geneset_tests.tsv: expected one row for {mark_label}/{gs}"
    return float(row.iloc[0][which])

# ---- Fig. 5A: significant-peak counts, counted from the DESeq2 result tables
fig, axes = plt.subplots(1, 4, figsize=(12, 3.2), sharey=True)
for ax, (mk, name) in zip(axes, MARKS):
    vals = []
    for x in CONTRASTS_5A:
        f = need(os.path.join(DESQ, mk, f"cells_Adult_vs_{x}.results.tsv"))
        r = pd.read_csv(f, sep="\t")
        vals.append(int((r["padj"] < 0.05).sum()))
    ax.bar(["Pooled", "LAP", "Dual", "Single"], vals, color="#4C72B0")
    ax.set_title(f"{name}\nAdult vs fetal, FDR<0.05", fontsize=9)
    ax.tick_params(labelsize=8)
axes[0].set_ylabel("Significant peaks")
save(fig, "fig5A_sig_peak_counts")

# ---- Fig. 5B: promoter boxplots (cells) ------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(12, 3.4), sharey=True)
for ax, (mk, name) in zip(axes, MARKS):
    d = ret[mk]
    prom = d[d.annotation.str.startswith("Promoter", na=False)]
    cu = prom[prom.deg_class == "colony_up"].cells_log2FC.dropna()
    hu = prom[prom.deg_class == "host_up"].cells_log2FC.dropna()
    ax.boxplot([cu, hu], labels=[f"colony-up\n(n={len(cu)})", f"host-up\n(n={len(hu)})"],
               showfliers=False)
    ax.axhline(0, color="grey", lw=0.8)
    p = pipe_p(name, "colony_up", "cells_p")
    ax.set_title(f"{name}: promoter peaks, cells\np = {p:.1e}", fontsize=9)
    ax.tick_params(labelsize=8)
axes[0].set_ylabel("log2FC (fetal vs adult)")
save(fig, "fig5B_promoter_boxplots")

# ---- Fig. 6A: retention scatters -------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(12, 3.2), sharey=True)
for ax, (mk, name) in zip(axes, MARKS):
    d = ret[mk].dropna(subset=["cells_log2FC", "tissue_log2FC"])
    if len(d) > SCATTER_MAX_POINTS:
        d = d.sample(SCATTER_MAX_POINTS, random_state=1)
    ax.scatter(d.cells_log2FC, d.tissue_log2FC, s=1, alpha=0.25)
    ax.axhline(0, color="grey", lw=0.8); ax.axvline(0, color="grey", lw=0.8)
    ax.set_title(f"{name}: retention (all peaks)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_xlabel("cells log2FC", fontsize=8)
axes[0].set_ylabel("tissue log2FC (colony vs host)")
save(fig, "fig6A_retention_scatters")

# ---- Concept Fig. A: promoter medians, cells vs tissue, pipeline p-values --
fig, axes = plt.subplots(2, 4, figsize=(12, 6.4), sharey="row",
                         gridspec_kw={"hspace": 0.75})
for j, (mk, name) in enumerate(MARKS):
    d = ret[mk]
    prom = d[d.annotation.str.startswith("Promoter", na=False)]
    for i, (col, lab) in enumerate([("cells_log2FC", "cells"), ("tissue_log2FC", "tissue")]):
        ax = axes[i][j]
        meds = [prom[prom.deg_class == c][col].median() for c in ["colony_up", "host_up"]]
        ax.bar(["colony-up", "host-up"], meds, color=["#0072B2", "#D55E00"])
        ax.axhline(0, color="grey", lw=0.8)
        p = pipe_p(name, "colony_up", "cells_p" if i == 0 else "tissue_p")
        ax.set_title(f"{name} {lab}\nmedian log2FC, p = {p:.1e}", fontsize=9)
        ax.tick_params(labelsize=8)
axes[0][0].set_ylabel("cells median log2FC"); axes[1][0].set_ylabel("tissue median log2FC")
fig.tight_layout()
save(fig, "figA_promoter_medians")

# ---- Concept Fig. B: H3K36me1 sensitivity sign-flip --------------------------
nd = pd.read_csv(need(os.path.join(SENS, "3m1", "tissue_nodedup.results.tsv")), sep="\t")
m = ret["3m1"].merge(nd[["chr", "start", "end", "log2FC_mle"]].rename(columns={"log2FC_mle": "tissue_nd"}),
                     on=["chr", "start", "end"], how="inner")
prom = m[m.annotation.str.startswith("Promoter", na=False)]
fig, ax = plt.subplots(figsize=(5.2, 3.6))
classes = ["colony_up", "host_up"]
dedup = [prom[prom.deg_class == c].tissue_log2FC.median() for c in classes]
nodup = [prom[prom.deg_class == c].tissue_nd.median() for c in classes]
x = np.arange(2); w = 0.36
ax.bar(x - w/2, dedup, w, label="deduplicated", color="#0072B2")
ax.bar(x + w/2, nodup, w, label="duplicates retained", color="#D55E00")
ax.axhline(0, color="grey", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(["colony-up promoters", "host-up promoters"], fontsize=8)
ax.set_ylabel("median tissue log2FC (colony vs host)")
ax.set_title("H3K36me1 tissue promoters: counting choice flips the sign\n(dedup p = 0.20; duplicates retained p = 0.002, anti-retention)", fontsize=9)
ax.legend(fontsize=8)
save(fig, "figB_h3k36me1_signflip", tight=True)

# ---- Concept Fig. C: QC, cells vs tissue ------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
for ax, col, lab, logy in [(axes[0], "frip", "FRiP", False),
                           (axes[1], "dup_pct", "Duplication %", False),
                           (axes[2], "unique_frags", "Unique fragments", True)]:
    for t, color in [("cells", "#0072B2"), ("tissue", "#D55E00")]:
        s = metrics[metrics.type == t]
        ax.scatter(s.mark, s[col], c=color, label=t, alpha=0.75)
    if logy: ax.set_yscale("log")
    ax.set_title(lab, fontsize=9); ax.tick_params(labelsize=8)
axes[0].legend(fontsize=8)
save(fig, "figC_qc_cells_vs_tissue")

# ---- Concept Fig. D: concordance summary -------------------------------------
sens_sum = pd.read_csv(need(os.path.join(SENS, "concordance_summary.tsv")), sep="\t")
rows = []
for mk, name in MARKS:
    d = ret[mk].dropna(subset=["cells_log2FC", "tissue_log2FC"])
    agree = float((np.sign(d.cells_log2FC) == np.sign(d.tissue_log2FC)).mean())
    rho = float(d[["cells_log2FC", "tissue_log2FC"]].corr(method="spearman").iloc[0, 1])
    srow = sens_sum[sens_sum.mark == mk].iloc[0]
    rows.append((name, agree, rho, float(srow.sign_agreement), float(srow.spearman_rho)))
fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
x = np.arange(4); w = 0.36
axes[0].bar(x - w/2, [r[1] for r in rows], w, label="deduplicated", color="#0072B2")
axes[0].bar(x + w/2, [r[3] for r in rows], w, label="duplicates retained", color="#D55E00")
axes[0].axhline(0.5, color="grey", lw=0.9, ls="--")
axes[0].set_xticks(x); axes[0].set_xticklabels([r[0] for r in rows], fontsize=8)
axes[0].set_ylabel("sign agreement"); axes[0].set_title("Genome-wide sign agreement (0.5 = chance)", fontsize=9)
axes[0].legend(fontsize=8)
axes[1].bar(x - w/2, [r[2] for r in rows], w, color="#0072B2")
axes[1].bar(x + w/2, [r[4] for r in rows], w, color="#D55E00")
axes[1].axhline(0, color="grey", lw=0.9)
axes[1].set_xticks(x); axes[1].set_xticklabels([r[0] for r in rows], fontsize=8)
axes[1].set_ylabel("Spearman rho"); axes[1].set_title("Cells vs tissue effect correlation", fontsize=9)
save(fig, "figD_concordance_summary")

print("ALL PANELS DONE ->", OUT_FIGDIR)
