#!/usr/bin/env python3

import json
import subprocess
from pathlib import Path
import sys


def run(cmd, shell=False):
    print(f"\n[RUNNING]\n{cmd if isinstance(cmd, str) else ' '.join(cmd)}\n")
    subprocess.run(cmd, shell=shell, check=True)


def require_file(path, label):
    path = Path(path)
    if not path.exists():
        sys.exit(f"ERROR: {label} not found: {path}")


def main(config_path):
    with open(config_path) as f:
        cfg = json.load(f)

    ref = Path(cfg["ref_fasta"])
    r1 = Path(cfg["r1"])
    r2 = Path(cfg["r2"])
    barcode_bed = Path(cfg["barcode_bed"])

    require_file(ref, "Reference FASTA")
    require_file(r1, "R1 FASTQ")
    require_file(r2, "R2 FASTQ")
    require_file(barcode_bed, "Barcode BED file")

    prefix = cfg.get("prefix", "sample")
    threads = int(cfg.get("threads", 8))
    min_qual = int(cfg.get("min_qual", 30))
    min_depth = int(cfg.get("min_depth", 10))

    outdir = Path(cfg.get("output_dir", "custom_reference_output"))
    outdir.mkdir(parents=True, exist_ok=True)

    trimmed_r1 = outdir / f"{prefix}.trim.R1.fastq.gz"
    trimmed_r2 = outdir / f"{prefix}.trim.R2.fastq.gz"
    bam = outdir / f"{prefix}.PA14.bam"
    vcf = outdir / f"{prefix}.filtered.vcf.gz"
    consensus = outdir / f"{prefix}.PA14_like_consensus.fasta"
    masked = outdir / f"{prefix}.PA14_like_consensus.barcode_masked.fasta"

    bakta_out = outdir / f"bakta_{prefix}"
    final_gbk = outdir / f"{prefix}.gbk"

    # 1. Trim reads
    run([
        "fastp",
        "-i", str(r1),
        "-I", str(r2),
        "-o", str(trimmed_r1),
        "-O", str(trimmed_r2),
        "--detect_adapter_for_pe",
        "--thread", str(threads),
        "--html", str(outdir / f"{prefix}.fastp.html"),
        "--json", str(outdir / f"{prefix}.fastp.json")
    ])

    # 2. Index reference
    run(["bwa", "index", str(ref)])

    # 3. Map reads to unsorted BAM
    unsorted_bam = outdir / f"{prefix}.PA14.unsorted.bam"

    map_cmd = (
        f"bwa mem -t {threads} {ref} {trimmed_r1} {trimmed_r2} | "
        f"samtools view -@ {threads} -bS - > {unsorted_bam}"
    )
    run(map_cmd, shell=True)

    # Sort BAM using older samtools-compatible syntax
    sort_prefix = outdir / f"{prefix}.PA14.sorted"
    run([
        "samtools", "sort",
        "-@", str(threads),
        str(unsorted_bam),
        str(sort_prefix)
    ])

    # Rename sorted BAM to expected name
    sorted_bam = outdir / f"{prefix}.PA14.sorted.bam"
    run(["mv", str(sorted_bam), str(bam)])

    run(["samtools", "index", str(bam)])

    # 4. Call and filter variants
    call_cmd = (
        f"bcftools mpileup -Ou -f {ref} {bam} | "
        f"bcftools call -mv -Ou | "
        f"bcftools filter -Oz "
        f"-i 'QUAL>={min_qual} && DP>={min_depth}' "
        f"-o {vcf}"
    )
    run(call_cmd, shell=True)

    run(["bcftools", "index", str(vcf)])

    # 5. Make PA14-like consensus
    consensus_cmd = (
        f"cat {ref} | bcftools consensus {vcf} > {consensus}"
    )
    run(consensus_cmd, shell=True)

    # 6. Mask barcode region
    run([
        "bedtools", "maskfasta",
        "-fi", str(consensus),
        "-bed", str(barcode_bed),
        "-fo", str(masked)
    ])

    # 7. Annotate with Bakta
    run([
        "bakta",
        "--force",
        "--db", cfg["bakta_db"],
        "--genus", cfg.get("genus", "Pseudomonas"),
        "--species", cfg.get("species", "aeruginosa"),
        "--strain", cfg.get("strain", prefix),
        "--output", str(bakta_out),
        "--prefix", prefix,
        str(masked)
    ])

    # 8. Copy Bakta GBFF to .gbk for Snippy
    gbff = bakta_out / f"{prefix}.gbff"
    if gbff.exists():
        run(["cp", str(gbff), str(final_gbk)])
        print(f"\nDONE. Use this as your Snippy reference:\n{final_gbk}\n")
    else:
        print("\nBakta finished, but I could not find the expected .gbff file.")
        print(f"Check this folder: {bakta_out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python make_custom_reference.py config.json")
    main(sys.argv[1])