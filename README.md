# vLoD - Variant Limit of Detection Tool

A high-performance Rust implementation of the vLoD tool for statistically assessing the detectability of alleles from variant call files (VCF) using matched sequencing data.

## Introduction

vLoD calculates the likelihood of observing each variant in the context of a given sequencing error rate, true positive rate, and false positive rate. This allows users to assign a detectability score to each variant and classify variants as detectable or non-detectable.

## Features

- **Detectability Scoring**: Computes a log-odds ratio for each variant representing the likelihood of it being a true positive
- **Parallel Processing**: Utilizes all available CPU cores for fast variant processing
- **VCF Integration**: Directly annotates input VCF files with detectability status and scores
- **Flexible Workflow**: Use the combined tool for simplicity, or separate tools for advanced pipelines
- **Compressed File Support**: Handles both plain and gzip-compressed VCF files

## Installation

### From Source

```bash
git clone https://github.com/akkusalper/vLoD.git
cd vLoD
cargo build --release
```

The compiled binaries will be available in `target/release/`.

### Requirements

- Rust 1.70+ (2021 edition)
- System libraries: zlib, libbz2, liblzma, htslib

On Ubuntu/Debian:
```bash
apt-get install zlib1g-dev libbz2-dev liblzma-dev libhts-dev
```

## Usage

### Combined Workflow (Recommended)

The `vlod` tool performs detectability analysis and VCF annotation in a single step:

```bash
vlod --input-vcf variants.vcf --input-bam alignments.bam --output annotated.vcf
```

This reads variants from the VCF, analyzes the BAM file to calculate detectability scores, and outputs an annotated VCF with two new INFO fields:
- `DET`: Detectability status (`Yes` or `No`)
- `DETS`: Detectability score (float)

### Advanced Two-Step Workflow

For pipelines requiring intermediate results, use the separate tools:

**Step 1: Analyze detectability**
```bash
lod_edit --input-vcf variants.vcf --input-bam alignments.bam --output results.tsv
```

**Step 2: Merge results into VCF**
```bash
merge_vcf_lod variants.vcf results.tsv annotated.vcf
```

### Command-Line Options

#### vlod / lod_edit

| Option | Default | Description |
|--------|---------|-------------|
| `--input-vcf` | Required | Path to input VCF file |
| `--input-bam` | Required | Path to input BAM file (index required) |
| `--output` | Required | Path to output file |
| `--TP` | 0.999 | Probability of true positive |
| `--FP` | 0.001 | Probability of false positive |
| `--SE` | 0.0001 | Probability of sequencing error |
| `--num-processes` | Auto | Number of parallel processes |
| `-v, --verbose` | Off | Enable verbose logging |
| `-d, --debug` | Off | Enable debug logging |
| `-f, --force` | Off | Overwrite output if exists |

#### merge_vcf_lod

```bash
merge_vcf_lod [OPTIONS] <VCF_FILE> <DETECTABILITY_FILE> <OUTPUT_FILE>
```

### BAM Index Requirement

The BAM file must have an accompanying index file. The tool automatically looks for:
- `<filename>.bam.bai`
- `<filename>.bai`

## Output Format

### TSV Output (lod_edit)

| Column | Description |
|--------|-------------|
| Chrom | Chromosome |
| Pos | Position (1-based) |
| Ref | Reference allele |
| Alt | Alternative allele |
| Detectability_Score | Log-odds ratio score |
| Detectability_Condition | "Detectable" or "Non-detectable" |
| Coverage | Total read depth at position |
| Variant_Reads | Number of reads supporting the variant |

### Annotated VCF Output

The output VCF includes two new INFO fields:
```
##INFO=<ID=DET,Number=1,Type=String,Description="Detectability status (Yes if detectable, No if non-detectable)">
##INFO=<ID=DETS,Number=1,Type=Float,Description="Detectability Score">
```

Example variant line:
```
chr1	100	.	A	T	.	PASS	DP=30;DET=Yes;DETS=3.5
```

## Detectability Classification

Variants are classified based on their detectability score:
- **Detectable**: Score ≥ 2.50
- **Non-detectable**: Score < 2.50

## Docker

A Docker image is available for containerized execution:

```bash
docker pull alperakkus/vlod:latest
```

See [Docker Hub](https://hub.docker.com/r/alperakkus/vlod) for the latest releases.

## Examples

Basic analysis with default parameters:
```bash
vlod --input-vcf sample.vcf --input-bam sample.bam --output sample_annotated.vcf
```

With custom probability parameters:
```bash
vlod --input-vcf sample.vcf --input-bam sample.bam --output sample_annotated.vcf \
    --TP 0.995 --FP 0.005 --SE 0.001
```

Verbose output with 8 processes:
```bash
vlod --input-vcf sample.vcf --input-bam sample.bam --output sample_annotated.vcf \
    --num-processes 8 --verbose
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use vLoD in your research, please cite the original work by Alper Akkus.
