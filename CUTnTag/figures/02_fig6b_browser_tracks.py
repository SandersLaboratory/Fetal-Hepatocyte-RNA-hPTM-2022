#!/usr/bin/env python3
"""
02_fig6b_browser_tracks.py — manuscript Fig. 6B: CPM signal tracks at the
candidate loci (default: Cacna1b for H3K27me1, Cacna1e for H3K27me3),
adult vs fetal cell populations, from the pipeline's CPM bigWigs.

For each locus: one panel with four population-mean tracks (Adult, LAP,
Dual, Single; mean of the two biological replicates, 100 bp bins) over the
gene body +/- FLANK_BP, plus a gene-model track from the GTF exons.
Output: fig6b_tracks_<gene>.pdf/.png in OUT_FIGDIR.

REVISED Oct 22, 2026 — three fixes over the original version:
  1. All four population tracks within a locus now share ONE y-axis
     scale (the maximum signal across the four populations at that
     locus, plus headroom), so track heights are directly comparable.
     Previously each track auto-scaled, which made cross-population
     comparisons misleading.
  2. Titles use the manuscript mark names (H3K27me1 / H3K27me3)
     instead of the pipeline mark codes (2m1 / 2m3).
  3. The x-axis label reads "Chr N position (strand ...)".

Requires pyBigWig (NOT in the base cutntag env). Preflight below checks and
prints the fix:
  conda activate cutntag && pip install pyBigWig
(or run inside any env that has pyBigWig + matplotlib + pandas.)

Gene coordinates are looked up from the GTF at runtime — nothing hardcoded.
"""
import os, sys, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import pyBigWig
except ImportError:
    sys.exit("pyBigWig is not installed in this environment.\n"
             "Fix:  conda activate cutntag && pip install pyBigWig")

# ============================ PARAMETERS — CONFIRM BEFORE RUNNING ============
CNT_OUT    = os.environ.get("CNT_OUT", "~/BrownRIH/cnt_output")  # <-- CONFIRM
BIGWIG_DIR = os.path.join(CNT_OUT, "bigwig")                                       # files: <sample>.bw
GTF        = os.environ.get("CNT_GTF", "~/BrownRIH/genome/Rattus_norvegicus.mRatBN7.2.110.gtf")  # env CNT_GTF overrides
OUT_FIGDIR = os.path.join(CNT_OUT, "figures_final")
# (gene symbol, mark code) — the two loci the Results text cites
LOCI = [("Cacna1b", "2m1"), ("Cacna1e", "2m3")]
MARK_NAMES = {"2m1": "H3K27me1", "2m3": "H3K27me3", "3m1": "H3K36me1", "3m2": "H3K36me2"}
# sample-id prefixes per population (samples are <prefix><rep>-<mark>)
POPULATIONS = [("Adult", ["A1", "A2"]), ("LAP", ["L1", "L2"]),
               ("Dual", ["D1", "D2"]), ("Single", ["S1", "S2"])]
FLANK_BP = 10_000
BIN_BP = 100
Y_HEADROOM = 1.05        # shared y-limit = locus max signal * this factor
COLORS = {"Adult": "#4C72B0", "LAP": "#55A868", "Dual": "#C44E52", "Single": "#8172B2"}
# =============================================================================

os.makedirs(OUT_FIGDIR, exist_ok=True)

def gene_locus(gtf_path, symbol):
    """Return (chrom, start, end, strand, [(exon_start, exon_end), ...]) for a gene."""
    pat = re.compile(r'gene_name "%s"' % re.escape(symbol))
    chrom = strand = None; gs = ge = None; exons = []
    with open(gtf_path) as fh:
        for line in fh:
            if not pat.search(line):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] == "gene":
                chrom, gs, ge, strand = f[0], int(f[3]), int(f[4]), f[6]
            elif f[2] == "exon":
                exons.append((int(f[3]), int(f[4])))
    if chrom is None:
        sys.exit(f"gene {symbol} not found in {gtf_path} (check gene_name spelling)")
    return chrom, gs, ge, strand, sorted(exons)

for symbol, mk in LOCI:
    chrom, gs, ge, strand, exons = gene_locus(GTF, symbol)
    lo, hi = max(1, gs - FLANK_BP), ge + FLANK_BP
    bins = np.arange(lo, hi, BIN_BP)
    centers = bins[:-1] + BIN_BP / 2

    # Pass 1: compute every population's mean signal so the four tracks
    # can share a single y-scale for this locus.
    signals = []
    for pop, prefixes in POPULATIONS:
        sig = np.zeros(len(centers)); n = 0
        for pre in prefixes:
            bw_path = os.path.join(BIGWIG_DIR, f"{pre}-{mk}.bw")
            if not os.path.exists(bw_path):
                sys.exit(f"MISSING BIGWIG: {bw_path}")
            bw = pyBigWig.open(bw_path)
            vals = np.array(bw.stats(chrom, lo, hi, type="mean", nBins=len(centers)), dtype=float)
            bw.close()
            sig += np.nan_to_num(vals); n += 1
        signals.append(sig / max(n, 1))
    ymax = max(float(s.max()) for s in signals) * Y_HEADROOM
    print(f"{symbol} ({MARK_NAMES[mk]}): shared y-limit 0-{ymax:.3f} CPM")

    # Pass 2: plot with the locked scale.
    fig, axes = plt.subplots(len(POPULATIONS) + 1, 1, figsize=(11, 5.2), sharex=True,
                             gridspec_kw={"height_ratios": [1] * len(POPULATIONS) + [0.35], "hspace": 0.12})
    for ax, (pop, _prefixes), sig in zip(axes, POPULATIONS, signals):
        ax.fill_between(centers, sig, color=COLORS[pop], alpha=0.75, lw=0)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(pop, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.set_xlim(lo, hi)
    axes[0].set_title(f"{symbol} — {MARK_NAMES[mk]} (CPM, population mean of 2 replicates)", fontsize=10)
    gax = axes[-1]
    for (s, e) in exons:
        gax.fill_between([s, e], [0, 0], [1, 1], color="black")
    gax.plot([gs, ge], [0.5, 0.5], color="black", lw=1)
    gax.set_yticks([]); gax.set_ylabel("gene", fontsize=9)
    gax.set_xlabel(f"Chr {chrom} position (strand {strand})", fontsize=9)
    for name in ("fig6b_tracks_" + symbol.lower(),):
        fig.savefig(os.path.join(OUT_FIGDIR, name + ".pdf"))
        fig.savefig(os.path.join(OUT_FIGDIR, name + ".png"), dpi=200)
    plt.close(fig)
    print("wrote fig6b_tracks_" + symbol.lower())
print("FIG 6B DONE ->", OUT_FIGDIR)
