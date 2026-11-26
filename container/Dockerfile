# BALLAD (B-ALL Allelic imbalance and Dosage): B-ALL copy-number calling from hg38 bulk RNA-seq BAMs.
# Everything except the hg38 reference FASTA is inside: tools, R packages, the B-ALL expression reference, the 1kGP
# panel site lists and the exon cache (resources/, ~290 MB).
#
# Build (from the repository root, so resources/ and R/ are in the context):
#   docker build --platform linux/amd64 -f container/Dockerfile -t ballad:0.1.2 .
# Run:
#   docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
#     -v /path/to/bams:/data:ro -v /path/to/hg38_dir:/ref:ro -v "$PWD/out":/out \
#     ballad:0.1.2 run --bam /data/S.bam --ref /ref/hg38.fa --outdir /out --threads 8

FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="BALLAD" \
      org.opencontainers.image.description="B-ALL karyotype, LOH and tumor fraction from hg38 bulk RNA-seq (allele counts at 1kGP SNPs + expression)" \
      org.opencontainers.image.version="0.1.2"

COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm -f /tmp/environment.yml

USER root
COPY ballad.py /opt/ballad/ballad.py
COPY R/baf_lib.R R/call_arms.R R/count_genes.R R/plot_genome.R R/run_numbat.R R/build_reference.R /opt/ballad/R/
COPY resources/ /opt/ballad/resources/
# all.sites.vcf.gz is the per-chromosome site lists joined (kept out of git: > 100 MB)
RUN cd /opt/ballad/resources/sites && \
    /opt/conda/bin/bcftools concat --no-version -Oz -o all.sites.vcf.gz $(for c in $(seq 1 22); do echo chr$c.sites.vcf.gz; done) && \
    /opt/conda/bin/bcftools index -f -t all.sites.vcf.gz
COPY container/config.json /opt/ballad/config.json
COPY container/entrypoint.sh /usr/local/bin/ballad
RUN chmod 755 /usr/local/bin/ballad /opt/ballad/ballad.py && \
    chmod -R a+rX /opt/ballad && \
    mkdir -p /ref /data /out && chmod 777 /out

# conda environment on PATH without activation, so plain `singularity exec` / `--user` runs work too
ENV PATH=/opt/conda/bin:$PATH \
    MPLCONFIGDIR=/tmp \
    HOME=/tmp \
    OMP_NUM_THREADS=1
WORKDIR /out
USER $MAMBA_USER
ENTRYPOINT ["/usr/local/bin/ballad"]
CMD ["--help"]
