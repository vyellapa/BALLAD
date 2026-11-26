#!/usr/bin/env bash
# Experimental: GLIMPSE2 genotype imputation + phasing of an RNA BAM against the lp_loh 1kGP reference (already split
# into GLIMPSE2 binary chunks by lp_loh/scripts/glimpse_prep_reference.sh). BAM mode: GLIMPSE2 computes genotype
# likelihoods from the reads. Output per chromosome: <outdir>/glimpse/<sample>.chr<N>.bcf with GT (phased), GP, DS.
# Chunks run concurrently, 1 thread each (single-sample GLIMPSE2_phase uses ~1 core).
# Usage: glimpse_rna.sh <bam> <sample> <sample outdir> [jobs] [chroms, default "1..22"]
set -euo pipefail
BAM=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); S=$2; SD=$3; JOBS=${4:-8}; CHROMS=${5:-$(seq 1 22)}
REF=${GLIMPSE_REF:?set GLIMPSE_REF to the GLIMPSE2 reference panel directory}
IMG=${GLIMPSE_IMG:-simrub/glimpse:v2.0.0-27-g0919952_20221207}
OUT=$(mkdir -p "$SD/glimpse/chunks" && cd "$SD/glimpse" && pwd)
BAMDIR=$(dirname "$BAM")
export BAM S REF IMG OUT BAMDIR

g() { docker run --rm --platform linux/amd64 -v "$REF:/ref:ro" -v "$BAMDIR:/bam:ro" -v "$OUT:/out" "$IMG" "$@"; }
export -f g

phase_chunk() {  # chrom, chunk index, input region
  local O=chunks/chr$1.$2.bcf BIN=1kGP.chr$1_${3//[:\-]/_}.bin
  [ -s "$OUT/$O" ] && [ -s "$OUT/$O.csi" ] && return 0
  g GLIMPSE2_phase --bam-file "/bam/$(basename "$BAM")" --ind-name "$S" --reference "/ref/$BIN" --threads 1 \
    --output "/out/$O" --log "/out/$O.log" > /dev/null
  [ -s "$OUT/$O" ] || { echo "GLIMPSE2_phase failed: $O (see $OUT/$O.log)" >&2; return 1; }
}
export -f phase_chunk

for C in $CHROMS; do cut -f1,3 "$REF/chunks.chr$C.txt" | sed "s/^/$C\t/"; done \
  | xargs -P "$JOBS" -L 1 bash -c 'phase_chunk "$0" "$1" "$2"'

for C in $CHROMS; do
  cut -f1 "$REF/chunks.chr$C.txt" | sed "s|^|/out/chunks/chr$C.|; s|$|.bcf|" > "$OUT/chunks/chr$C.list"
  g GLIMPSE2_ligate --input "/out/chunks/chr$C.list" --output "/out/$S.chr$C.bcf" --threads 2 > "$OUT/ligate.chr$C.log" 2>&1
  [ -s "$OUT/$S.chr$C.bcf" ] || { echo "GLIMPSE2_ligate failed: chr$C (see $OUT/ligate.chr$C.log)" >&2; exit 1; }
  [ -s "$OUT/$S.chr$C.bcf.csi" ] || bcftools index -f "$OUT/$S.chr$C.bcf"
done
echo "glimpse done: $OUT"
