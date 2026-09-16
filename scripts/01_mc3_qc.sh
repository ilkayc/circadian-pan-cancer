#!/bin/bash

set -euo pipefail

DATA="data/mc3.v0.2.8.PUBLIC.xena.gz"
OUT="results/mc3_qc"

mkdir -p "$OUT"

echo "=== MC3 QC START ==="
date

echo ""
echo "1. Total mutation events:"
gzip -dc "$DATA" | tail -n +2 | wc -l

echo ""
echo "2. Unique samples:"
gzip -dc "$DATA" | tail -n +2 | cut -f1 | sort -u | wc -l

echo ""
echo "3. Unique genes:"
gzip -dc "$DATA" | tail -n +2 | cut -f7 | sort -u | wc -l

echo ""
echo "4. Effect classes:"
gzip -dc "$DATA" | tail -n +2 | cut -f8 | sort | uniq -c | sort -nr > "$OUT/effect_counts.txt"
cat "$OUT/effect_counts.txt"

echo ""
echo "5. Gene mutation counts:"
gzip -dc "$DATA" | tail -n +2 | cut -f7 | sort | uniq -c | sort -nr > "$OUT/gene_mutation_counts.txt"
head -30 "$OUT/gene_mutation_counts.txt"

echo ""
echo "6. Sample mutation counts:"
gzip -dc "$DATA" | tail -n +2 | cut -f1 | sort | uniq -c | sort -nr > "$OUT/sample_mutation_counts.txt"
head -20 "$OUT/sample_mutation_counts.txt"

echo ""
echo "7. Chromosome distribution:"
gzip -dc "$DATA" | tail -n +2 | cut -f2 | sort | uniq -c | sort -k1,1nr > "$OUT/chromosome_counts.txt"
cat "$OUT/chromosome_counts.txt"

echo ""
echo "=== MC3 QC COMPLETE ==="
date
