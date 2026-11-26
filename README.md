# BALLAD

**B-ALL Allelic imbalance and Dosage**: chromosome-level copy-number calling for B-ALL from bulk RNA-seq: BAM → allele counts at known 1kGP SNPs (no variant
calling) → genotype-free allelic imbalance + expression against a diploid B-ALL reference (numbat's expression model),
segmented jointly → karyotype, ploidy group, LOH and tumor fraction. No phasing, no variant calling.

Input: an hg38 RNA-seq BAM (or, for archived samples, a HaplotypeCaller VCF + HTSeq counts). The Docker image holds
everything except the hg38 FASTA. Validation sample: **HH-1**, a high-hyperdiploid diagnostic bone-marrow RNA-seq
sample (clinical subtype hyperdiploid by cytogenetics / FISH / CMA); its data are not included.

## Requirements
- `bcftools`, `samtools`, `bgzip`, `tabix`, `Rscript`, Python ≥ 3.9 (standard library only)
- R packages: numbat (+ hahmmr), data.table, dplyr, Matrix, ggplot2, glue, optparse, org.Hs.eg.db, AnnotationDbi
  (`prepare` only), GenomicAlignments, Rsamtools, rtracklayer (`count` step only)
- cellsnp-lite (default pileup): on PATH, or the numbat container `pkharchenkolab/numbat-rbase` (~21 GB compressed, amd64)
  via docker or a Singularity `.sif` (`cellsnp_container` in `config.json`). `--pileup bcftools` needs no container.
- `resources/sites/all.sites.vcf.gz` is built from the per-chromosome site lists on first run (too large for git)
- Paths to the hg38 FASTA, Ensembl 102 GTF (gene set of the RNAseqCNV reference), 1kGP hg38 panel (site list and allele frequencies) and
  RNAseqCNV source are in `config.json`

## Docker

Everything except the hg38 FASTA is in the image (`container/`, details in `container/README.md`):

```bash
docker build --platform linux/amd64 -f container/Dockerfile -t ballad:0.1.0 .
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
  -v /path/to/bams:/data:ro -v /path/to/hg38_dir:/ref:ro -v "$PWD/out":/out \
  ballad:0.1.0 run --bam /data/S.bam --ref /ref/hg38.fa --outdir /out --threads 8
```

## Use
```bash
# per sample: BAM (hg38, indexed) + the reference it was aligned to; everything else comes from config.json
./ballad.py run --bam S.bam --ref hg38.fa --outdir out --threads 8
#   sample name = BAM file name (or --sample); gene counts are made from the BAM (or --counts S.HTSeq);
#   the BAM and FASTA are checked to be hg38; the one-time resources (below) are built on first use

# no BAM (archived samples): allele counts from the HaplotypeCaller VCF (variant sites only) + HTSeq counts; same
# caller, flagged vcf_input; LOH where the reference allele is kept can be missed or shortened
./ballad.py run --vcf S.HaplotypeCaller.vcf --counts S.HTSeq --outdir out

# one-time resources, explicitly (e.g. with in-house diploid samples as the expression reference)
./ballad.py prepare --threads 8 [--standard-counts in_house_diploid_htseq_files.txt]
```
`--steps count,pileup,numbat` reruns a subset. The FASTA is only read by `--pileup bcftools` and the build
check; cellsnp-lite counts at the panel sites without it.

### On a Linux server, inside the numbat image
The numbat image (`pkharchenkolab/numbat-rbase:v1.4.2`, or a `.sif` built from it) already has Python 3, R + numbat,
cellsnp-lite, samtools, bgzip and tabix, so the whole run can happen inside it. The image lacks bcftools,
org.Hs.eg.db and GenomicAlignments, so:
- run `prepare` on a machine that has them and copy `resources/` along with the code (`resources/ens_to_symbol.tsv`
  replaces org.Hs.eg.db at run time);
- give `--counts S.HTSeq` (the clinical HTSeq file); the `count` step needs GenomicAlignments.

Copy the 1kGP panel (`panel_dir`, 5.1 GB) and set its path in `config.server.json`. Then:
```bash
singularity exec -B /path/to/bams,/path/to/htseq,/path/to/panel_hg38,/path/to/out numbat.sif \
  python3 /path/to/ballad/ballad.py --config /path/to/ballad/config.server.json \
  run --bam S.bam --counts S.HTSeq --outdir /path/to/out --threads 8
```
Tested with the docker image (same calls as the native run on HH-1).

## Outputs (`out/<sample>/`)
| file | content |
|---|---|
| `calls_summary.tsv` | **main result (BAM input)**: chromosome number, ploidy group, gains / losses / whole-chromosome LOH, LOH regions, focal events (likely allele-specific expression or expression artefacts; not counted), blast fraction with its 95% segment-bootstrap interval (`blast_fraction_lo` / `_hi`), model parameters, alternative-ploidy score gap, flags |
| `calls_arms.tsv` | per arm: allelic imbalance `theta`, baseline-corrected log2 expression, call, fraction of the arm it covers, events, numbat's call |
| `calls_segments.tsv` | segments: `theta` (with support interval) and log2 expression, observed vs expected for the called state (total copies `c`, minor copies `minor`); both-allele rate observed vs expected (`both_ratio`, `both_ratio_exp`); the blast fraction this segment implies on its own (`p_implied`, 1-SE `_lo` / `_hi`, from BAF or, for 2:2 gains, expression) |
| `<sample>_blast_fraction.png` | per-segment implied blast fraction against the joint fit and bootstrap band (clonal events line up), and the fit score over p |
| `blast_fraction_profile.tsv`, `blast_fraction_boot.tsv` | best score at each p (baseline profiled out); each bootstrap resample's p and the fitted segments it left out |
| `summary.tsv` | numbat's own expression-only summary (cross-check, not the calls): estimated chromosome number, ploidy group, gains / losses, flags. Expression alone cannot see CN-LOH or tell 2:2 from 3:1; on HH-1 it reads 48 (−19, −22 that the allele data rejects) |
| `karyotype.tsv` | per autosome: numbat's expression-only call, fraction of the chromosome it covers, sub-chromosomal events, expression ratio `phi`, both-allele rate |
| `<sample>_genome.png` | two panels: expression log2 ratio (faint = single genes, coloured = 25-gene bins) and BAF (faint = single sites with ≥ 10 reads, coloured = genotype-free bins at 0.5 ± theta), coloured by the `calls_segments.tsv` call, segments as lines, called events labelled; blast fraction and interval in the subtitle |
| `gene_counts.tsv` | gene counts from the BAM (HTSeq format; absent when `--counts` is given) |
| `allele_counts.tsv.gz` | panel sites with reads: REF/ALT counts, 1kGP AF |
| `rnaseqcnv_snv.tsv` | the same counts as RNAseqCNV's custom SNV table (chr start depth maf) |
| `numbat/` | numbat expression-only bulk table (gene logFC; read by the caller and the plot) and segments |

## Method notes
- **Pileup**: cellsnp-lite bulk mode at the panel sites, as numbat's `pileup_and_phase.R --bulk` (MAPQ ≥ 20; duplicates
  skipped unless `keep_duplicates` in `config.json`); alternative `bcftools mpileup` (MAPQ ≥ 20, BQ ≥ 20, no BAQ). Counting at fixed known sites (not variant calling) keeps sites where one allele
  is lost and makes every sample comparable.
- **No phasing** (Eagle2 removed Nov 2026): the caller's allele model never used phase, RNA phase holds only within
  ~50 kb, and the phased BAF track only added noise to the plot. numbat now runs expression-only; it still needs
  allele rows, so hets by its genotype rule are passed with a placeholder phase (identical segments and logFC to the
  former phased run on HH-1). Removing it drops the Eagle / docker dependency on macOS.
- **Expression reference**: RNAseqCNV's 44 diploid B-ALL samples (or in-house via `--standard-counts`); numbat fits each
  sample as a mixture of them. MHC is masked by numbat.
- **B-ALL prior** (`call_arms.R --ball_prior`): the balanced 4:2 gain (tetrasomy 21, the hallmark of high
  hyperdiploidy) is charged as one event, like a single-copy gain.
- **Chromosome number** = 46 + whole-chromosome changes (gain +1, balanced gain +2, loss −1, CN-LOH 0). Autosomes only;
  sex chromosomes assumed 2 - check X with RNAseqCNV.
- **Flags**: `possible_doubled_hypodiploid` (≥ 5 whole-chromosome CN-LOH: a doubled hypodiploid looks hyperdiploid or
  diploid by expression and shows up as widespread LOH), `low_allele_data` (< 5,000 het sites),
  `partial_chromosome_events` (a chromosome's main call covers < 80% of it).
- **Allele + expression caller** (`R/call_arms.R`; model in `R/baf_lib.R`; with `--vcf` input it sees variant sites
  only and is flagged `vcf_input`):
  - *Allelic imbalance without genotypes or phase*: every panel SNP with reads, at any depth, is a mixture of
    hom-ref / het / hom-alt with the 1kGP prior, het prior scaled by a fitted sample-wide het rate `h`; het reads
    0.5 Bin(d, 0.5+θ) + 0.5 Bin(d, 0.5−θ). Calling hets first drops low-depth hets that show one allele, which hides
    LOH (numbat's allele HMM misses HH-1's 16p CN-LOH for this reason). Bins of ~20 expected hets (0.3-1.5 Mb).
  - *Overdispersion*: likelihood curves are tempered by κ (median per-bin LR statistic vs the chi-square(1) median;
    4.7 on HH-1). Balanced RNA regions sit at θ0 ≈ 0.05, not 0 (allele-specific expression); expected θ is
    √(θ² + θ0²).
  - *Segmentation*: optimal partitioning per chromosome on θ curves + gene log2 expression jointly, ≥ 3 bins per segment.
  - *States*: (total, minor) copies in {(1,0) loss, (2,0) CN-LOH, (2,1), (3,1) gain, (3,0), (4,2) balanced double
    gain, (4,1)}, blast fraction and expression baseline fitted jointly. Expression cannot fix absolute ploidy, so the
    best solution with a different mean copy number is reported (`alt_mean_cn`, `alt_score_gap`) and flagged as
    `ploidy_ambiguous` when close (near-haploid vs doubled hypodiploid).
  - *Segment score* (adopted from HAPLOID `joint_fit.py`): the segment's own θ log-likelihood curve read at the predicted
    θ (keeps the curve's shape; sd floored at `--theta_sd` 0.02) + expression z² + a both-allele term (sites with each
    allele on ≥ max(2, 10%) of reads, observed vs expected from 1kGP AF, het rate and depth, normalised to balanced
    segments; it reuses the same reads, so it is a second view of LOH, not new data) + `--state_penalty` per copy
    changed (PURPLE event distance: 3:1 and 1:0 = 1, 2:0 and 4:1 = 2, 3:0 = 3; `--ball_prior` charges the balanced 4:2
    gain as 1). With one flat charge per state, a 3:1 gain at p = 1 and a 4:1 at p ≈ 0.5 cost the same and p rested on
    the 16p CN-LOH alone (`--penalty flat` restores it).
  - *Blast fraction*: grid p 0.05-1 (step 0.01) × expression baseline, fitted on segments ≥ 5 bins and ≥ `--fit_min_mb`
    5 Mb (a 3.1 Mb ASE segment at 11q12-13 otherwise pulled bootstrap resamples to p ≈ 0.48). 95% interval from 200
    segment-bootstrap resamples (`--boot`): how much p depends on which events drive it. It is a share of RNA, not of
    cells; blasts and normal marrow cells may make different amounts of RNA. A CN-LOH that is germline homozygosity
    says nothing about p - check whether its `p_implied` agrees with the gains.
  - *Focal*: non-neutral segments < 5 bins or < 3 Mb are listed (`focal`), not counted - in RNA these are mostly
    mono-allelic expression (11p15.5 imprinted cluster, 5q31 protocadherins on HH-1). Balanced double gains rest
    on expression alone and need ≥ 10 bins; shorter ones are listed as focal expression changes (1q21 pericentromeric
    genes with DRAGEN-BAM counts against the STAR-based reference).
  - *LOH*: CN-LOH can be germline homozygosity (runs of homozygosity) as well as somatic; confirm with remission DNA.
  - Tested and not adopted (scripts in `experimental/`): GLIMPSE2 genotype prior and GLIMPSE phasing (no gain; RNA phase
    concordance 88-92% < 10 kb, ~55% at 200 kb-1 Mb), beta-binomial overdispersion (not identifiable), binned phased BAF (Eagle2; removed).
- **Both-allele rate** (pileup input only): observed vs expected sites showing both alleles, from 1kGP AF and depth.
  Low values mean hets have dropped out (LOH at high blast %).

## Validation status
Tested end-to-end on one high-hyperdiploid sample (HH-1, VCF input): +6, +14, +17, +21×2, 51 chromosomes.
Pileup + genotype checked for mechanics on a DNA BAM (7.6 GB WGS) with both tools: cellsnp-lite (numbat image, amd64
emulation, 8 threads) 13 min, 4.92 M sites, 531 k hets; bcftools 4.91 M sites, 516 k hets; same genotype class at 99.2% of
shared sites. cellsnp-lite counts ~1.5 more reads where they differ (no base-quality filter, overlapping mates counted
twice). Full BAM run (cellsnp-lite) on HH-1 RNA BAM (15 GB): 7.6 min, 364 k sites with reads, 27 k hets phased by Eagle2 at the time (24.6 k from
the VCF; phasing since removed); same karyotype as the VCF run; both-allele rate 0.75-1.26 on every autosome. The allele + expression caller
on the same BAM: +6, +14, +17, +21×2 (balanced), 51 chromosomes, blast fraction ~1.0, plus 16p CN-LOH (0-14 Mb, missed by
numbat; somatic vs germline unresolved) and 5 focal allele-specific-expression regions. The caller has not yet seen a loss,
a hypodiploid or a doubled hypodiploid. BAM-only run (`run --bam --ref`, gene counts from the DRAGEN BAM, reverse
strand detected; 8 min 51 s, 8 threads): identical calls. The clinical HTSeq files come from a different (STAR-style)
alignment, so BAM-derived counts differ gene by gene (log correlation 0.84); the expression reference should ideally be
counted the same way as the samples (`prepare --standard-counts`). Without phasing (Nov 2026): calls and blast
fraction unchanged; numbat's expression-only cross-check reads 48. VCF input (variant sites only, 49.5 k sites) through
the caller: same 51-chromosome karyotype, 16p CN-LOH 0.1-12.7 Mb, blast fraction 0.97 (0.88-1.00). Blast fraction on HH-1: 0.97 (95% bootstrap
0.89-1.00), calls unchanged. In-silico dilutions of the same sample with its own normal cells
(`validation/sim_admix.R`, results in `validation/dilution results (not included: patient-derived)`, 3 seeds each; expected p = 0.97 × the
dilution): recovered 0.77 at 0.78, 0.48-0.49 at 0.48 and 0.30-0.31 at 0.29, every interval containing the expected
value; at 0.19 the estimate rises to 0.27-0.38 and the 16p CN-LOH is miscalled as a loss. Event calls fail before p
does: all events at 0.8 (one seed lost +6 to a 6q piece read as 2:2), +17 lost in 1/3 at 0.5 and in 3/3 at 0.3 (+6 and +14 also in 1/3 at 0.3),
+21×2 kept at 0.3 only with `--ball_prior`, nothing called at 0.15. The dilution model assumes normal cells express at
the reference level and share the tumor's germline genotype. Before relying on the pipeline: run it on RNA BAMs with known
karyotypes (diploid, high hyperdiploid, low hypodiploid / near-haploid, doubled hypodiploid) and calibrate the B-ALL prior (`--ball_prior`),
the CN-LOH flag threshold and the both-allele rate.

## License and third-party components

BALLAD is MIT-licensed (`LICENSE`). Please also cite what it builds on:

| used for | software / data | license |
|---|---|---|
| B-ALL diploid expression reference, gene weights (`resources/ball_reference.rds`) | RNAseqCNV (Bařinka J, et al. Leukemia 2022) | MIT |
| expression model (bulk mode) | numbat, hahmmr (Gao T, et al. Nat Biotechnol 2023) | MIT |
| allele counting | cellsnp-lite (Huang X, Huang Y. Bioinformatics 2021); samtools / bcftools / htslib | MIT / MIT, GPL |
| SNP sites and allele frequencies (`resources/sites/`) | 1000 Genomes Project phase 3, hg38 | open |
| purity bootstrap, event-distance prior, clonal tumor-fraction check | adapted from HAPLOID (PURPLE event distance) | MIT |
