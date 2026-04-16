#!/bin/bash
set -euo pipefail

THREADS=${1:-47}

for bam in *.sorted.bam; do
    base="${bam%.sorted.bam}"
    out="${base}.called.bam"

    if [[ -s "$out" ]]; then
        echo "Skipping existing $out"
        continue
    fi

    echo "Calling mods for $bam -> $out"
    modkit call-mods -t "$THREADS" "$bam" "$out"
    samtools index -@ "$THREADS" "$out"
done
