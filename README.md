# make_custom_reference

`make_custom_reference` builds a sample-specific, PA14-like reference genome from paired-end Illumina reads. It trims reads, maps them to a supplied reference FASTA, calls high-confidence variants, applies those variants to make a consensus FASTA, masks a configured barcode region, annotates the masked consensus with Bakta, and copies the resulting GenBank file for downstream tools such as Snippy.

## Repository contents

| File | Purpose |
| --- | --- |
| `make_custom_reference.py` | Main pipeline script. Reads a JSON config, runs the external bioinformatics tools, and writes all outputs to the configured output directory. |
| `config.json` | Example configuration file with input paths, filtering thresholds, Bakta metadata, and output settings. |
| `barcode_region.bed` | Example BED file describing the region to mask in the consensus genome. |

## Workflow

The script performs these steps:

1. Validate that the reference FASTA, paired FASTQ files, and barcode BED file exist.
2. Trim paired-end reads with `fastp`.
3. Index the reference with `bwa index`.
4. Align reads to the reference with `bwa mem` and convert the alignment to BAM with `samtools view`.
5. Sort and index the BAM with `samtools`.
6. Call variants with `bcftools mpileup` and `bcftools call`.
7. Filter variants by configurable `QUAL` and `DP` thresholds.
8. Create a consensus FASTA with `bcftools consensus`.
9. Mask the barcode region with `bedtools maskfasta`.
10. Annotate the masked consensus with `bakta`.
11. Copy Bakta's `.gbff` output to a final `.gbk` file suitable for use as a Snippy reference.

## Requirements

Install the following command-line tools and make sure they are available on your `PATH`:

- Python 3
- `fastp`
- `bwa`
- `samtools`
- `bcftools`
- `bedtools`
- `bakta`

You also need a local Bakta database and must provide its path in `config.json` with the `bakta_db` setting.

## Configuration

Copy or edit `config.json` before running the pipeline. The important settings are:

| Setting | Description |
| --- | --- |
| `ref_fasta` | Path to the reference genome FASTA. |
| `r1` / `r2` | Paths to paired-end FASTQ files. |
| `prefix` | Sample/output prefix used for generated filenames. Defaults to `sample` if omitted. |
| `threads` | Number of threads for supported tools. Defaults to `8` if omitted. |
| `barcode_bed` | BED file containing the barcode region(s) to mask. |
| `min_qual` | Minimum variant quality used by `bcftools filter`. Defaults to `30` if omitted. |
| `min_depth` | Minimum variant depth used by `bcftools filter`. Defaults to `10` if omitted. |
| `bakta_db` | Path to the Bakta database directory. |
| `genus`, `species`, `strain` | Metadata passed to Bakta. |
| `output_dir` | Directory for all generated outputs. Defaults to `custom_reference_output` if omitted. |

BED coordinates are zero-based and half-open, as expected by `bedtools maskfasta`.

## How to run

From the repository root, run:

```bash
python3 make_custom_reference.py config.json
```

For a new sample, create a new config file and pass it to the script:

```bash
cp config.json my_sample.config.json
# Edit paths, prefix, barcode BED, Bakta metadata, and output_dir as needed.
python3 make_custom_reference.py my_sample.config.json
```

## Outputs

If `prefix` is `my_barcode_pool` and `output_dir` is `PA14wa_BC_output`, the pipeline creates outputs such as:

| Output | Description |
| --- | --- |
| `PA14wa_BC_output/my_barcode_pool.trim.R1.fastq.gz` | Trimmed R1 reads. |
| `PA14wa_BC_output/my_barcode_pool.trim.R2.fastq.gz` | Trimmed R2 reads. |
| `PA14wa_BC_output/my_barcode_pool.PA14.bam` | Sorted alignment BAM. |
| `PA14wa_BC_output/my_barcode_pool.PA14.bam.bai` | BAM index. |
| `PA14wa_BC_output/my_barcode_pool.filtered.vcf.gz` | Filtered variant calls. |
| `PA14wa_BC_output/my_barcode_pool.PA14_like_consensus.fasta` | Consensus FASTA with sample variants applied. |
| `PA14wa_BC_output/my_barcode_pool.PA14_like_consensus.barcode_masked.fasta` | Consensus FASTA with barcode region masked. |
| `PA14wa_BC_output/bakta_my_barcode_pool/` | Bakta annotation output directory. |
| `PA14wa_BC_output/my_barcode_pool.gbk` | Final GenBank file intended for use as a Snippy reference. |

## Notes and limitations

- The script runs external commands directly and exits if any command fails.
- The reference is indexed in place by `bwa index`, so the directory containing `ref_fasta` must be writable.
- The bundled `config.json` contains environment-specific absolute paths; update them before running on another system.
- The script currently uses an older `samtools sort` output-prefix style. If your installed `samtools` requires `-o`, update that command or use a compatible version.
