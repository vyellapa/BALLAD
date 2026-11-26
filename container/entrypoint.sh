#!/usr/bin/env bash
# Entry point of the ballad image: runs the pipeline with the config baked into the image.
#   docker run ... ballad run --bam /data/S.bam --ref /ref/hg38.fa --outdir /out --threads 8
#   docker run ... ballad run --vcf /data/S.HaplotypeCaller.vcf --counts /data/S.HTSeq --outdir /out
#   docker run ... ballad check        (tool versions and built-in resources)
set -euo pipefail
APP=/opt/ballad
case "${1:-}" in
  ""|-h|--help|help)
    exec python3 "$APP/ballad.py" --config "$APP/config.json" --help ;;
  check)
    echo "ballad image"
    python3 --version
    R --version | head -1
    Rscript -e 'for (p in c("numbat", "hahmmr", "GenomicAlignments", "Rsamtools", "rtracklayer", "data.table", "ggplot2", "patchwork")) cat(sprintf("%-18s %s\n", p, as.character(packageVersion(p))))'
    echo "cellsnp-lite $(cellsnp-lite --version 2>&1 | head -1)"
    samtools --version | head -1
    bcftools --version | head -1
    echo "resources: $(ls "$APP/resources" | tr '\n' ' ')"
    echo "panel sites: $(zcat "$APP/resources/sites/all.sites.vcf.gz" | grep -vc '^#')"
    exit 0 ;;
  bash|sh)
    exec "$@" ;;
  *)
    exec python3 "$APP/ballad.py" --config "$APP/config.json" "$@" ;;
esac
