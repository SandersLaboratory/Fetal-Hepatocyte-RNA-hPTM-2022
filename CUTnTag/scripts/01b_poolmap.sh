#!/bin/bash
# 01b_poolmap.sh — reconstruct sequencing pool/lane mapping from FASTQ headers.
# Azenta never provided the pool map for the index-conflict split (2 pools of 16
# for the 32 cell samples; tissue ran as its own project). Illumina headers
# (@INST:RUN:FLOWCELL:LANE:TILE:X:Y) recover it exactly.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUT="$OUTDIR/qc/pool_map.tsv"
echo -e "sample\tinstrument\trun\tflowcell\tlane" > "$OUT"
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample r1 r2 rest; do
  # NOTE: do NOT use zcat ... | head -1 here. Under `set -o pipefail`, head
  # closing the pipe gives zcat SIGPIPE (exit 141), and `set -e` then kills
  # the whole script on the first sample. Python reads one line with no pipe.
  H=$(python3 -c 'import gzip,sys; sys.stdout.write(gzip.open(sys.argv[1],"rt").readline().strip())' \
      "$FASTQ_DIR/$r1")
  echo -e "${sample}\t$(echo "$H" | cut -d: -f1 | sed 's/^@//')\t$(echo "$H" | cut -d: -f2)\t$(echo "$H" | cut -d: -f3)\t$(echo "$H" | cut -d: -f4)" >> "$OUT"
done

echo "=== samples per flowcell:lane (each = one sequencing pool) ==="
tail -n +2 "$OUT" | awk -F'\t' '{k=$4":lane"$5; a[k]=a[k]" "$1} END{for(k in a) print k" ->"a[k]}' | sort
echo "(pools inferred from headers override samples.tsv in script 04)"
