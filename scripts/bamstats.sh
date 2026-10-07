#!/bin/bash
# usage: bamstats.sh BAM BED THREADS OUT_TXT OUT_SAMTOOLS_STATS
set -euo pipefail
BAM=$1; BED=$2; T=$3; OUT=$4; ST=$5
samtools view -@ "$T" -F 260 "$BAM" | awk '{q=$5; n++; s+=q; if(q>=20)a++; if(q>=30)b++}
  END{printf "n_primary_mapped=%d\nmapq_mean=%.2f\nmapq20=%d\nmapq30=%d\n", n, s/n, a, b}' > "$OUT"
echo "on_target_reads=$(samtools view -@ "$T" -c -F 260 -L "$BED" "$BAM")" >> "$OUT"
samtools stats -@ "$T" "$BAM" > "$ST"
