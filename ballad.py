#!/usr/bin/env python3
"""B-ALL copy-number calling from bulk RNA-seq: BAM -> allele counts at 1kGP SNPs -> genotype-free allelic imbalance +
expression (numbat's B-ALL expression model) -> copy-number calls, karyotype and blast fraction. No phasing.

  prepare : once per install. B-ALL expression reference (R/build_reference.R) and the SNP site lists: 1kGP hg38 biallelic
            SNPs (MAF >= 1%) inside numbat's gene bodies, chr1-22 (resources/sites/).
  run     : per sample. Needs only the BAM (hg38, aligned to the same build as --ref / config fasta); the sample name
            defaults to the BAM file name, and prepare runs first if its resources are missing.
    count    gene counts from the BAM (R/count_genes.R, htseq-count union mode, Ensembl 102 = the gene set of the
             RNAseqCNV reference, strandedness detected) unless --counts gives an HTSeq file. -> gene_counts.tsv
    pileup   allele counts at the panel sites. Backends: cellsnp-lite (default; on PATH, else inside the numbat container,
             docker or Singularity - config.json cellsnp_container) or bcftools mpileup (--pileup bcftools; MAPQ/BQ
             filters, no BAQ). Duplicates skipped unless keep_duplicates.
             --vcf instead of --bam: counts from a VCF's FORMAT/AD (HaplotypeCaller fallback when the BAM is gone;
             variant sites only, so LOH where the reference allele is kept can be missed or shortened).
             -> allele_counts.tsv.gz: CHROM POS REF ALT AF REF_N ALT_N DP OTH (DP = REF_N + ALT_N), panel alleles.
    numbat   R/run_numbat.R: expression-only numbat with the B-ALL reference (gene logFC for the caller, and numbat's
             own expression-only calls as a cross-check: karyotype.tsv, summary.tsv); then R/call_arms.R (allele +
             expression calls per segment / arm / chromosome, blast fraction: calls_*.tsv) and R/plot_genome.R
             (<sample>_genome.png). Allelic imbalance comes from every panel site with reads (R/baf_lib.R), so no
             genotyping or phasing step is needed (Eagle2 phasing was dropped: RNA phase holds only within ~50 kb and
             the calls never used it).
  Also written: rnaseqcnv_snv.tsv (chr start depth maf) - RNAseqCNV's custom SNV table, so it can use the same counts.

Usage:
  ballad.py run --bam S.bam --ref hg38.fa --outdir out [--threads 8]
  ballad.py run --bam S.bam --outdir out [--sample S] [--counts S.HTSeq] [--pileup bcftools] [--steps ...]
  ballad.py run --vcf S.HaplotypeCaller.vcf --counts S.HTSeq --outdir out     (no BAM; variant sites only)
  ballad.py prepare [--threads 8] [--standard-counts list.txt]               (runs automatically if needed)
"""
import argparse
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHROMS = [str(c) for c in range(1, 23)]
STEPS = ["count", "pileup", "numbat"]
HG38_CHR1 = 248956422


def log(msg):
    print(f"[ballad] {msg}", file=sys.stderr, flush=True)


def load_config(path):
    cfg = json.loads(Path(path).read_text())
    for k in ("fasta", "gtf", "panel_dir", "rnaseqcnv_dir", "resources_dir"):
        p = Path(os.path.expanduser(cfg[k]))
        cfg[k] = p if p.is_absolute() else HERE / p
    cfg["sites_dir"] = cfg["resources_dir"] / "sites"
    return cfg


def sh(cmd, log_file=None):
    """Run a bash pipeline with pipefail; stderr to log_file if given."""
    err = open(log_file, "a") if log_file else None
    try:
        subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, stderr=err)
    finally:
        if err:
            err.close()


def need(*tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        sys.exit(f"missing tools on PATH: {', '.join(missing)}")


def pmap(fn, items, threads):
    with ThreadPoolExecutor(max_workers=threads) as ex:
        return list(ex.map(fn, items))


# ---------------------------------------------------------------- prepare

def prepare(cfg, args):
    need("Rscript", "bcftools", "bgzip", "tabix")
    res, sites = cfg["resources_dir"], cfg["sites_dir"]
    sites.mkdir(parents=True, exist_ok=True)
    cmd = (f"Rscript {HERE / 'R' / 'build_reference.R'} --rnaseqcnv_dir {cfg['rnaseqcnv_dir']} --out_dir {res}"
           + (f" --standard_counts {args.standard_counts}" if args.standard_counts else ""))
    log("building B-ALL expression reference")
    sh(cmd)

    # gene bodies, merged, as 0-based BED
    genes = {}
    with open(res / "genes.tsv") as f:
        next(f)
        for line in f:
            c, s, e, _ = line.rstrip("\n").split("\t")
            genes.setdefault(c, []).append((int(s) - 1, int(e)))
    with open(res / "genes_merged.bed", "w") as out:
        for c in [f"chr{x}" for x in CHROMS]:
            merged = []
            for s, e in sorted(genes.get(c, [])):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            out.writelines(f"{c}\t{s}\t{e}\n" for s, e in merged)

    def one(chrom):
        panel = cfg["panel_dir"] / cfg["panel_pattern"].format(chrom=chrom)
        vcf, tsv, tgt = (sites / f"chr{chrom}.{x}" for x in ("sites.vcf.gz", "sites.tsv.gz", "targets.tsv.gz"))
        sh(f"bcftools view --no-version -G -T {res / 'genes_merged.bed'} -Oz -o {vcf} {panel} && bcftools index -f -t {vcf}")
        sh(f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' {vcf} | bgzip > {tsv}")
        sh(f"bgzip -dc {tsv} | cut -f1,2 | bgzip > {tgt} && tabix -f -s1 -b2 -e2 {tgt}")
        n = int(subprocess.run(f"bgzip -dc {tsv} | wc -l", shell=True, capture_output=True, text=True).stdout)
        log(f"chr{chrom}: {n} panel SNPs in genes")
        return n

    log("extracting panel SNPs inside genes")
    total = sum(pmap(one, CHROMS, args.threads))
    all_vcf = sites / "all.sites.vcf.gz"      # single site list for cellsnp-lite -R
    sh(f"bcftools concat --no-version -Oz -o {all_vcf} " + " ".join(str(sites / f'chr{c}.sites.vcf.gz') for c in CHROMS)
       + f" && bcftools index -f -t {all_vcf}")
    log(f"done: {total} sites; resources in {res}")


# ---------------------------------------------------------------- run: pileup

def load_sites(cfg, chrom):
    sites = {}
    with gzip.open(cfg["sites_dir"] / f"chr{chrom}.sites.tsv.gz", "rt") as f:
        for line in f:
            _, pos, ref, alt, af = line.rstrip("\n").split("\t")
            sites[int(pos)] = (ref, alt, af)
    return sites


def bam_contig_style(bam):
    hdr = subprocess.run(["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True).stdout
    names = {l.split("SN:")[1].split("\t")[0] for l in hdr.splitlines() if l.startswith("@SQ")}
    if "chr1" in names:
        return "chr"
    if "1" in names:
        return ""
    sys.exit(f"{bam}: no chr1 / 1 contig in the header")


def pileup_bcftools(cfg, args, sd):
    need("bcftools", "samtools")
    pc = cfg["pileup"]
    prefix = bam_contig_style(args.bam)
    fai = Path(str(cfg["fasta"]) + ".fai")
    if fai.exists() and f"{prefix}1" not in {l.split("\t")[0] for l in fai.read_text().splitlines()}:
        sys.exit(f"{cfg['fasta']} contig names do not match the BAM ({prefix or 'no '}chr prefix)")
    skip = "UNMAP,SECONDARY,QCFAIL" if pc["keep_duplicates"] else "UNMAP,SECONDARY,QCFAIL,DUP"
    tmp = sd / "tmp"; tmp.mkdir(exist_ok=True)

    def one(chrom):
        targets = cfg["sites_dir"] / f"chr{chrom}.targets.tsv.gz"
        if not prefix:                                   # BAM without chr prefix: strip it from the targets
            plain = tmp / f"chr{chrom}.targets.tsv"
            sh(f"bgzip -dc {targets} | sed 's/^chr//' > {plain}")
            targets = plain
        cmd = (f"bcftools mpileup -f {shlex.quote(str(cfg['fasta']))} -r {prefix}{chrom} -T {targets} "
               f"-a FORMAT/AD -q {pc['min_mapq']} -Q {pc['min_baseq']} -B -I -d {pc['max_depth']} "
               f"--ns {skip} -Ou {shlex.quote(str(args.bam))} "
               f"| bcftools query -f '%POS\\t%REF\\t%ALT\\t[%AD]\\n'")
        out = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
        if out.returncode:
            raise RuntimeError(f"mpileup chr{chrom} failed:\n{out.stderr[-2000:]}")
        sites, rows = load_sites(cfg, chrom), []
        for line in out.stdout.splitlines():
            pos, ref, alts, ad = line.split("\t")
            site = sites.get(int(pos))
            if site is None or ref != site[0]:
                continue
            ad = [int(x) if x != "." else 0 for x in ad.split(",")]
            alts = alts.split(",")
            ref_n = ad[0]
            alt_n = ad[1 + alts.index(site[1])] if site[1] in alts else 0
            if ref_n + alt_n > 0:
                rows.append((f"chr{chrom}", int(pos), site[0], site[1], site[2], ref_n, alt_n,
                             ref_n + alt_n, sum(ad) - ref_n - alt_n))
        log(f"pileup chr{chrom}: {len(rows)} sites with reads")
        return rows

    return [r for rows in pmap(one, CHROMS, args.threads) for r in rows]


def cellsnp_prefix(cfg, mounts):
    """Command prefix that runs cellsnp-lite: on PATH if installed, else inside the numbat container (cellsnp_container in
    config.json: docker image, or a Singularity .sif on HPC). mounts: {host_dir: container_dir}."""
    cc = cfg["cellsnp_container"]
    if shutil.which(cc["cellsnp"]):
        return cc["cellsnp"], {h: h for h in mounts}
    if cc.get("sif"):
        need("singularity")
        binds = ",".join(f"{h}:{c}" for h, c in mounts.items())
        return f"singularity exec -B {binds} {cc['sif']} {cc['cellsnp']}", mounts
    need("docker")
    if subprocess.run(["docker", "image", "inspect", cc["image"]], capture_output=True).returncode:
        sys.exit(f"docker image {cc['image']} not found locally (~21 GB compressed); `docker pull {cc['image']}`, "
                 f"or set cellsnp_container.sif, or use --pileup bcftools")
    vols = " ".join(f"-v {h}:{c}" for h, c in mounts.items())
    return f"docker run --rm --platform linux/amd64 {vols} {cc['image']} {cc['cellsnp']}", mounts


def pileup_cellsnp(cfg, args, sd):
    # cellsnp-lite bulk mode (one BAM, no cell barcodes) at the panel sites given with -R, as numbat's
    # pileup_and_phase.R --bulk does. cellSNP.base.vcf.gz INFO: AD = alt reads, DP = ref + alt reads, OTH = other alleles.
    # cellsnp-lite's own -p does not speed this mode up (one bulk BAM with -R: its threads share about one core), so
    # each chromosome runs as its own cellsnp-lite process, --threads at a time, largest chromosomes first.
    pc = cfg["pileup"]
    out = sd / "cellsnp"; out.mkdir(exist_ok=True)
    bam = args.bam.resolve()
    prefix, m = cellsnp_prefix(cfg, {str(bam.parent): "/bam", str(cfg["sites_dir"]): "/sites", str(out): "/out"})
    flags = "UNMAP,SECONDARY,QCFAIL" if pc["keep_duplicates"] else "UNMAP,SECONDARY,QCFAIL,DUP"

    def one(chrom):
        (out / f"chr{chrom}").mkdir(exist_ok=True)
        sh(f"{prefix} -s {m[str(bam.parent)]}/{bam.name} -I {args.sample} -O {m[str(out)]}/chr{chrom} "
           f"-R {m[str(cfg['sites_dir'])]}/chr{chrom}.sites.vcf.gz -p 1 --minMAF 0 --minCOUNT 1 "
           f"--minMAPQ {pc['min_mapq']} --exclFLAG {flags} --UMItag None --cellTAG None --gzip",
           sd / "logs" / f"cellsnp.chr{chrom}.log")

    pmap(one, CHROMS, args.threads)
    afs = {c: load_sites(cfg, c) for c in CHROMS}
    rows = []
    for part in CHROMS:
        vcf = out / f"chr{part}" / "cellSNP.base.vcf.gz"
        if not vcf.exists():
            sys.exit(f"cellsnp-lite wrote no {vcf}; see {sd / 'logs' / f'cellsnp.chr{part}.log'}")
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c, pos, _, ref, alt, _, _, info = line.split("\t")[:8]
                chrom = c.removeprefix("chr")
                site = afs.get(chrom, {}).get(int(pos))
                if site is None or (ref, alt) != site[:2]:
                    continue
                kv = dict(x.split("=") for x in info.strip().split(";") if "=" in x)
                alt_n, dp, oth = int(kv["AD"]), int(kv["DP"]), int(kv["OTH"])
                if dp == 0:  # only non-panel alleles seen
                    continue
                rows.append((f"chr{chrom}", int(pos), ref, alt, site[2], dp - alt_n, alt_n, dp, oth))
    return rows


def counts_from_vcf(cfg, args, sd):
    need("bcftools")
    q = subprocess.run(["bash", "-o", "pipefail", "-c",
                        f"bcftools view -m2 -M2 -v snps {shlex.quote(str(args.vcf))} -Ou "
                        f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n'"],
                       capture_output=True, text=True, check=True).stdout
    afs = {c: load_sites(cfg, c) for c in CHROMS}
    rows = []
    for line in q.splitlines():
        c, pos, ref, alt, ad = line.split("\t")
        chrom = c.removeprefix("chr")
        site = afs.get(chrom, {}).get(int(pos))
        if site is None or (ref, alt) != site[:2] or "." in ad:
            continue
        ref_n, alt_n = (int(x) for x in ad.split(",")[:2])
        if ref_n + alt_n > 0:
            rows.append((f"chr{chrom}", int(pos), ref, alt, site[2], ref_n, alt_n, ref_n + alt_n, 0))
    return rows


def step_pileup(cfg, args, sd):
    if args.vcf:
        log("allele counts from VCF (fallback; sites with an allele lost are missing)")
        rows = counts_from_vcf(cfg, args, sd)
    elif args.pileup == "cellsnp":
        rows = pileup_cellsnp(cfg, args, sd)
    else:
        rows = pileup_bcftools(cfg, args, sd)
    rows.sort(key=lambda r: (int(r[0][3:]), r[1]))
    with gzip.open(sd / "allele_counts.tsv.gz", "wt") as f:
        f.write("CHROM\tPOS\tREF\tALT\tAF\tREF_N\tALT_N\tDP\tOTH\n")
        f.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    with open(sd / "rnaseqcnv_snv.tsv", "w") as f:           # RNAseqCNV custom table: sites with alt reads
        f.write("chr\tstart\tdepth\tmaf\n")
        f.writelines(f"{r[0]}\t{r[1]}\t{r[7]}\t{r[6] / r[7]:.4f}\n" for r in rows if r[6] > 0)
    log(f"allele_counts.tsv.gz: {len(rows)} panel sites with reads")


# ---------------------------------------------------------------- run: counts, numbat, calls

def step_count(cfg, args, sd):
    if args.user_counts:
        log(f"using supplied gene counts {args.counts}")
        return
    if not args.bam:
        sys.exit("--vcf input needs --counts (gene counts are made from a BAM)")
    need("Rscript")
    sh(f"Rscript {HERE / 'R' / 'count_genes.R'} --bam {shlex.quote(str(args.bam))} --gtf {cfg['gtf']} "
       f"--exons_rds {cfg['resources_dir'] / 'exons_ensembl102.rds'} --out {args.counts} --threads {args.threads}",
       sd / "logs" / "count.log")
    log((sd / "logs" / "count.log").read_text().strip().splitlines()[-1])


def check_build(cfg, args):
    # the SNP panel, sites and annotation are hg38: refuse a BAM (or --ref) from another build
    hdr = subprocess.run(["samtools", "view", "-H", str(args.bam)], capture_output=True, text=True, check=True).stdout
    sq = dict(re.findall(r"@SQ\tSN:(\S+)\tLN:(\d+)", hdr))
    chr1 = sq.get("chr1") or sq.get("1")
    if chr1 is None or int(chr1) != HG38_CHR1:
        sys.exit(f"{args.bam}: chr1 length {chr1}, expected {HG38_CHR1} (hg38/GRCh38)")
    fai = Path(str(cfg["fasta"]) + ".fai")
    if fai.exists():
        ref = dict(l.split("\t")[:2] for l in fai.read_text().splitlines())
        rc1 = ref.get("chr1") or ref.get("1")
        if rc1 is None or int(rc1) != HG38_CHR1:
            sys.exit(f"{cfg['fasta']}: chr1 length {rc1}, expected {HG38_CHR1} (hg38/GRCh38)")
    elif args.pileup == "bcftools":
        sys.exit(f"{cfg['fasta']}.fai missing (samtools faidx {cfg['fasta']})")


def step_numbat(cfg, args, sd):
    need("Rscript")
    nc = cfg["numbat"]
    sh(f"Rscript {HERE / 'R' / 'run_numbat.R'} --sample {args.sample} --counts {shlex.quote(str(args.counts))} "
       f"--sample_dir {sd} --reference {cfg['resources_dir'] / 'ball_reference.rds'} "
       f"--t {nc['t']} --cnloh_flag {nc['cnloh_flag']} "
       f"--allele_source {'vcf' if args.vcf else 'pileup'}",
       sd / "logs" / "numbat.log")
    # allele + expression caller; with --vcf it sees variant sites only (flagged in calls_summary.tsv)
    sh(f"Rscript {HERE / 'R' / 'call_arms.R'} --sample {args.sample} --sample_dir {sd} "
       f"--allele_source {'vcf' if args.vcf else 'pileup'}", sd / "logs" / "calls.log")
    sh(f"Rscript {HERE / 'R' / 'plot_genome.R'} --sample {args.sample} --sample_dir {sd}", sd / "logs" / "plot.log")


def ensure_all_sites(cfg):
    """resources/sites/all.sites.vcf.gz (cellsnp-lite's site list) is the per-chromosome site files joined; it is over
    GitHub's 100 MB limit, so it is not in the repository and is rebuilt here (and in the Docker build) when missing."""
    sites = cfg["sites_dir"]
    all_vcf = sites / "all.sites.vcf.gz"
    parts = [sites / f"chr{c}.sites.vcf.gz" for c in CHROMS]
    if all_vcf.exists() or not all(x.exists() for x in parts):
        return
    need("bcftools")
    log(f"building {all_vcf} from the per-chromosome site lists (once)")
    sh(f"bcftools concat --no-version -Oz -o {all_vcf} " + " ".join(str(x) for x in parts) + f" && bcftools index -f -t {all_vcf}")


def run(cfg, args):
    if not (args.bam or args.vcf):
        sys.exit("run needs --bam (or --vcf as a fallback)")
    if args.ref:
        cfg["fasta"] = args.ref.resolve()
    if args.bam:
        args.bam = args.bam.resolve()
        check_build(cfg, args)
    if not args.sample:
        args.sample = re.sub(r"\.(bam|cram|vcf|vcf\.gz|bcf)$", "", (args.bam or args.vcf).name)
    ensure_all_sites(cfg)
    if not (cfg["sites_dir"] / "all.sites.vcf.gz").exists() or not (cfg["resources_dir"] / "ball_reference.rds").exists():
        log("resources missing: running prepare first")
        prepare(cfg, argparse.Namespace(threads=args.threads, standard_counts=None))
    sd = Path(args.outdir).resolve() / args.sample
    (sd / "logs").mkdir(parents=True, exist_ok=True)
    args.user_counts = args.counts is not None
    args.counts = args.counts.resolve() if args.counts else sd / "gene_counts.tsv"
    steps = args.steps.split(",") if args.steps else STEPS
    for step in STEPS:
        if step in steps:
            log(f"== {step}")
            globals()[f"step_{step}"](cfg, args, sd)
    if "numbat" in steps:
        sys.stdout.write((sd / "summary.tsv").read_text())
        if (sd / "calls_summary.tsv").exists():
            sys.stdout.write((sd / "calls_summary.tsv").read_text())
    log(f"outputs in {sd}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(HERE / "config.json"))
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="build the B-ALL reference and SNP site lists (once)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--standard-counts", help="file listing in-house diploid HTSeq count files (default: RNAseqCNV's)")
    r = sub.add_parser("run", help="call CNVs for one sample")
    r.add_argument("--bam", type=Path, help="RNA-seq BAM (hg38), indexed")
    r.add_argument("--ref", type=Path, help="hg38 FASTA the BAM was aligned to (default: config.json fasta)")
    r.add_argument("--sample", help="sample name (default: BAM file name)")
    r.add_argument("--vcf", type=Path, help="fallback when no BAM: allele counts from this VCF's FORMAT/AD")
    r.add_argument("--counts", type=Path, help="HTSeq gene counts (Ensembl IDs); default: counted from the BAM")
    r.add_argument("--outdir", required=True)
    r.add_argument("--threads", type=int, default=4)
    r.add_argument("--pileup", choices=["cellsnp", "bcftools"], default="cellsnp")
    r.add_argument("--steps", help=f"comma-separated subset of {','.join(STEPS)} (default: all)")
    args = ap.parse_args()
    cfg = load_config(args.config)
    {"prepare": prepare, "run": run}[args.cmd](cfg, args)


if __name__ == "__main__":
    main()
