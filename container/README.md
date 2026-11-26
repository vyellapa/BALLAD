# ballad container

One image runs the whole analysis on hg38 B-ALL RNA-seq: gene counts, allele counts at 1kGP SNPs, numbat's expression
model, the allele + expression caller (karyotype, LOH, tumor fraction) and the plots. Everything is inside the image
except the hg38 reference FASTA:

| inside the image | |
|---|---|
| tools | Python 3.12, R 4.4.3, numbat 1.5.2, hahmmr 1.0.0, GenomicAlignments 1.42, Rsamtools 2.22, rtracklayer 1.66, cellsnp-lite 1.2.3, samtools / bcftools / htslib 1.21 |
| `resources/` | B-ALL diploid expression reference (RNAseqCNV's 44 samples), Ensembl 102 exon cache and gene symbols, 5.4 M 1kGP hg38 SNP sites inside genes (site lists only; the 5 GB panel is not needed at run time) |

Image: `ballad:0.1.1`, linux/amd64, 1.2 GB. `docker run ... ballad:0.1.1 check` prints the versions.

## Build

From the repository root (the build context must contain `ballad.py`, `R/` and `resources/`):

```bash
docker build --platform linux/amd64 -f container/Dockerfile -t ballad:0.1.1 .
```

The conda packages come from conda-forge and bioconda only (`container/environment.yml`, no defaults channel).

## Run

```bash
mkdir -p out
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
  -v /path/to/bams:/data:ro \
  -v /path/to/hg38_dir:/ref:ro \
  -v "$PWD/out":/out \
  ballad:0.1.1 run --bam /data/S.bam --ref /ref/hg38.fa --outdir /out --threads 8
```

- The BAM must be hg38 / GRCh38 and indexed (`.bai` next to it); other builds are refused.
- `--ref`: the FASTA the BAM was aligned to, with its `.fai`. The default pileup (cellsnp-lite) does not read it, so
  the run also works without it (the build is then checked from the BAM header only); `--pileup bcftools` needs it.
- `--counts /data/S.HTSeq` uses existing HTSeq gene counts (Ensembl IDs) instead of counting from the BAM.
- No BAM (archived samples): `run --vcf /data/S.HaplotypeCaller.vcf --counts /data/S.HTSeq --outdir /out`. Variant
  sites only; flagged `vcf_input` in `calls_summary.tsv`.
- `--sample NAME` sets the output folder name (default: the input file name).

Outputs go to `out/<sample>/`: `calls_summary.tsv` (main result: chromosome number, ploidy group, gains, losses, LOH,
tumor fraction with 95% interval, flags), `calls_segments.tsv`, `calls_arms.tsv`, `<sample>_genome.png`,
`<sample>_blast_fraction.png`, numbat's expression-only cross-check (`summary.tsv`, `karyotype.tsv`) and `logs/`.
See the main README for every column.

## HPC (Singularity / Apptainer)

Without a registry, move the image as a file:

```bash
docker save ballad:0.1.1 | gzip > ballad_0.1.1.tar.gz      # on the build machine
singularity build ballad_0.1.1.sif docker-archive://ballad_0.1.1.tar.gz   # on the cluster
singularity run --cleanenv --bind /path/to/bams:/data,/path/to/hg38_dir:/ref,$PWD/out:/out \
  ballad_0.1.1.sif run --bam /data/S.bam --ref /ref/hg38.fa --outdir /out --threads 8
```

## Tested (Nov 2026, Docker 29 on Apple silicon, amd64 emulation)

- BAM steps on a chr21-22 slice of an hg38 BAM: gene counting (strandedness detection), cellsnp-lite and bcftools
  pileups (117.0 k vs 116.9 k sites), with and without the FASTA mounted; an hg19 BAM is refused.
- Not yet run in the container on a full RNA-seq BAM end to end; do that on the first cluster sample.
