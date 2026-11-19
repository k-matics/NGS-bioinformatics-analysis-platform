#!/usr/bin/env bash
# Build a DRAGEN-style ALT-masked GRCh38 reference from the 1000 Genomes analysis set.
#
# What it does (nothing existing is modified or overwritten):
#   1. checks the source FASTA + mask BED and records their SHA-256
#   2. bedtools maskfasta: replaces the regions in the DRAGEN mask BED with N
#      -> a NEW FASTA; the source FASTA is kept unchanged
#   3. verifies contig names/lengths are identical to the source
#   4. samtools faidx, samtools dict (new MD5s of the masked sequence), bwa index
#   5. writes <fasta>.masking.json (provenance; the app re-verifies it)
#
# Usage: build_masked_reference.sh <source.fa> <mask.bed> <out.fa> [mask_bed_url] [source_url]
set -euo pipefail

SRC=${1:?source FASTA}; BED=${2:?mask BED}; OUT=${3:?output FASTA}
BED_URL=${4:-}; SRC_URL=${5:-}

for f in "$SRC" "$BED"; do [[ -s $f ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done
for f in "$OUT" "$OUT.fai" "$OUT.bwt" "$OUT.masking.json" "${OUT%.fa}.dict"; do
  [[ -e $f ]] && { echo "ERROR: $f already exists - refusing to overwrite" >&2; exit 1; }
done
for t in bedtools samtools bwa sha256sum; do command -v $t >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done

log() { echo "[$(date '+%F %T')] $*"; }
OUTDIR=$(dirname "$OUT")
BED_COPY="$OUTDIR/$(basename "$BED")"
[[ -e $BED_COPY ]] || cp "$BED" "$BED_COPY"

log "SHA-256 of inputs"
SRC_SHA=$(sha256sum "$SRC" | cut -d' ' -f1); BED_SHA=$(sha256sum "$BED" | cut -d' ' -f1)
echo "  source $SRC_SHA"; echo "  mask   $BED_SHA"

[[ -s $SRC.fai ]] || { log "samtools faidx (source)"; samtools faidx "$SRC"; }

log "Checking every mask BED contig exists in the source"
missing=$(cut -f1 "$BED" | sort -u | comm -23 - <(cut -f1 "$SRC.fai" | sort -u) | wc -l)
[[ $missing -eq 0 ]] || { echo "ERROR: $missing mask contigs not in the source FASTA" >&2; exit 1; }

log "bedtools maskfasta -> $OUT"
bedtools maskfasta -fi "$SRC" -bed "$BED" -fo "$OUT"
samtools faidx "$OUT"

log "Verifying contig names and lengths are unchanged"
diff <(cut -f1,2 "$SRC.fai") <(cut -f1,2 "$OUT.fai") >/dev/null || { echo "ERROR: contigs differ after masking" >&2; exit 1; }

log "Sequence dictionary (MD5 of the masked sequences)"
samtools dict -a GRCh38 -s "Homo sapiens" -u "file://$OUT" -o "${OUT%.fa}.dict" "$OUT"

MASKED_BP=$(awk '{s+=$3-$2} END{print s}' "$BED"); NIV=$(grep -cv '^#' "$BED")

log "bwa index (about 1 hour, ~5 GB RAM)"
bwa index -a bwtsw "$OUT" 2>&1 | tail -3

log "Writing provenance $OUT.masking.json"
cat > "$OUT.masking.json" <<EOF
{
  "method": "DRAGEN-style ALT masking (bedtools maskfasta with the DRAGMAP hg38_alt_mask.bed)",
  "source_fasta": "$SRC",
  "source_url": "$SRC_URL",
  "source_sha256": "$SRC_SHA",
  "mask_bed": "$BED",
  "mask_bed_url": "$BED_URL",
  "mask_bed_sha256": "$BED_SHA",
  "mask_bed_copy": "$BED_COPY",
  "intervals": $NIV,
  "masked_bp": $MASKED_BP,
  "bedtools": "$(bedtools --version)",
  "samtools": "$(samtools --version | head -1)",
  "bwa": "$(bwa 2>&1 | grep Version | cut -d' ' -f2)",
  "created": "$(date -Iseconds)"
}
EOF
log "Done: $OUT"
