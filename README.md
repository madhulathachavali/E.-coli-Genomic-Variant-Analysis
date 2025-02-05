# E.coli-Genomic-Variant-Analysis

### 
This project focuses on processing 15 E. coli samples to generate and analyze VCF files using Nextflow. The workflow uses command-line tools such as Nextflow, bcftools, and the SRA Toolkit (https://www.ncbi.nlm.nih.gov/sra?term=SAMN02604091)

## Sample Information 

| Sample ID       | Run Accession  | Instrument              | Strategy       | # of Spots | # of Bases | Size   | Description |
|----------------|---------------|-------------------------|---------------|------------|------------|--------|-------------|
| EK-bisulfite   | SRR13921543    | Illumina MiSeq         | Bisulfite-Seq | 2M         | 601.8M     | 297.5Mb | Bisulfite-seq for methylase specificity |
| EK02-dnaseq    | SRR13921545    | Illumina NovaSeq 6000  | WGS           | 5.1M       | 1G         | 308.2Mb | Whole Genome Sequencing |
| EK01-dnaseq    | SRR13921546    | Illumina NovaSeq 6000  | WGS           | 4.9M       | 970.4M     | 288.6Mb | Whole Genome Sequencing |
| EK02-rims-3h-rep | SRR13921549 | Illumina NovaSeq 6000  | RIMS-seq (3h NaOH) | 7.4M | 1.5G | 443.5Mb | RIMS-seq 3-hour incubation |
| EK02-rims-1h-rep | SRR13921550 | Illumina NovaSeq 6000  | RIMS-seq (1h NaOH) | 7.2M | 1.4G | 429Mb | RIMS-seq 1-hour incubation |
| EK02-rims-10min-rep | SRR13921551 | Illumina NovaSeq 6000 | RIMS-seq (10min NaOH) | 8.5M | 1.7G | 511.5Mb | RIMS-seq 10-minute incubation |
| EK02-control-3h-rep | SRR13921552 | Illumina NovaSeq 6000 | Control TE (3h) | 7.2M | 1.4G | 435.2Mb | Control TE 3-hour incubation |
| EK02-control-1h-rep | SRR13921554 | Illumina NovaSeq 6000 | Control TE (1h) | 7.2M | 1.4G | 434.5Mb | Control TE 1-hour incubation |
| EK02-control-10min-rep | SRR13921555 | Illumina NovaSeq 6000 | Control TE (10min) | 7.4M | 1.5G | 447.5Mb | Control TE 10-minute incubation |
| EK01-rims-3h | SRR13921556 | Illumina NovaSeq 6000 | RIMS-seq (3h NaOH) | 8.1M | 1.6G | 482.9Mb | RIMS-seq 3-hour incubation |
| EK01-rims-1h | SRR13921557 | Illumina NovaSeq 6000 | RIMS-seq (1h NaOH) | 6.4M | 1.3G | 384Mb | RIMS-seq 1-hour incubation |
| EK01-rims-10min | SRR13921558 | Illumina NovaSeq 6000 | RIMS-seq (10min NaOH) | 7.3M | 1.5G | 442.6Mb | RIMS-seq 10-minute incubation |
| EK01-control-3h | SRR13921559 | Illumina NovaSeq 6000 | Control TE (3h) | 5.8M | 1.2G | 347Mb | Control TE 3-hour incubation |
| EK01-control-1h | SRR13921560 | Illumina NovaSeq 6000 | Control TE (1h) | 7.1M | 1.4G | 430.2Mb | Control TE 1-hour incubation |
| EK01-control-10min | SRR13921561 | Illumina NovaSeq 6000 | Control TE (10min) | 6M | 1.2G | 358.3Mb | Control TE 10-minute incubation |
| EcoliK12-CCS-15kb | SRR10971019 | PacBio Sequel | WGS (HiFi Reads) | 95,514 | 1.4G | 1.1Gb | PacBio HiFi sequencing with 15kb size selection |



### Objectives
- Perform quality control (QC) on raw sequencing data.
- Convert raw data to SAM, BAM, and eventually VCF format.

### Steps
#### 1. Sample Data
- Samples are based on the E. coli K-12 MG1655 strain.
- Raw sequencing data is retrieved using the SRA Toolkit.
#### 2. Reference Genome
- A reference genome is indexed for alignment and variant calling processes.
#### 3. Quality Control and Preprocessing
- Quality control is performed using FastQC and other QC tools.
- Trimmomatic is used for trimming and preprocessing the reads.
#### 4. Data Processing
- Raw reads are aligned to the reference genome.
- Data is sequentially converted from SAM → BAM → VCF format.
#### 5. Analysis
- The generated VCF files are analyzed for Variant calling, Allele frequency, and Phylogenetic relationships among the samples
#### Tools and Technologies
- Nextflow: Workflow management for reproducible data processing pipelines.
- SRA Toolkit: Retrieval of raw sequencing data.
- FastQC & Trimmomatic: Quality control and read trimming.
- 
## VCF Analysis and Downstream Processing

## 📌 Prerequisites

- [GATK](https://gatk.broadinstitute.org)
- [SnpEff](http://snpeff.sourceforge.net)
- [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)
- [Samtools](http://www.htslib.org)
- [IGV](https://software.broadinstitute.org/software/igv/home)

## Step 1 Validate VCF File

```bash
gatk ValidateVariants -V /Users/madhuchavali/Desktop/Ecoli/variants.vcf -R /Users/madhuchavali/Desktop/Ecoli/sequence.fasta
```

## Step 2 Filter Variants (GATK VariantFiltration)

```bash
gatk VariantFiltration \  
  -V /Users/madhuchavali/Desktop/Ecoli/variants.vcf \  
  --filter-expression "QD < 2.0 || FS > 60.0" \  
  --filter-name "LowQuality" \  
  -O /Users/madhuchavali/Desktop/Ecoli/filtered_variants.vcf
```

## Step 3 Annotate Variants using SnpEff

```bash
java -jar /Users/madhuchavali/Desktop/Ecoli/snpEff/snpEff.jar NC_000913_3 \  
  /Users/madhuchavali/Desktop/Ecoli/filtered_variants.vcf > /Users/madhuchavali/Desktop/Ecoli/filtered_variants_annotated.vcf
```

## Step 4 Extract High-Impact Variants

```bash
grep "HIGH" /Users/madhuchavali/Desktop/Ecoli/filtered_variants_annotated.vcf > /Users/madhuchavali/Desktop/Ecoli/high_impact_variants.vcf
```

## Step 5 Annotate with dbSNP (SnpSift)

```bash
java -jar /Users/madhuchavali/Desktop/Ecoli/snpEff/SnpSift.jar annotate \  
  -v /Users/madhuchavali/Desktop/Ecoli/dbSNP.vcf.gz \  
  /Users/madhuchavali/Desktop/Ecoli/filtered_variants_annotated.vcf \  
  > /Users/madhuchavali/Desktop/Ecoli/filtered_variants_annotated_with_dbSNP.vcf
```

## Step 6 Prepare BAM Files for IGV

### Move BAM Files

```bash
mkdir -p /Users/madhuchavali/Desktop/Ecoli/Ecoli_BAM_Files
mv ~/Downloads/*.bam /Users/madhuchavali/Desktop/Ecoli/Ecoli_BAM_Files/
mv ~/Downloads/*.bai /Users/madhuchavali/Desktop/Ecoli/Ecoli_BAM_Files/
```

### Check BAM File Validity

```bash
samtools quickcheck -v /Users/madhuchavali/Desktop/Ecoli/Ecoli_BAM_Files/*.bam
```

### Index BAM Files

```bash
for bam in /Users/madhuchavali/Desktop/Ecoli/Ecoli_BAM_Files/*.sorted.bam; do
    samtools index "$bam"
done
```

### Confirm Indexing

```bash
ls -lh /Users/madhuchavali/Desktop/Ecoli/Ecoli_BAM_Files/ | grep ".bai"
```

## Step 7: Load BAM and VCF Files into IGV

## Output Files
| Step | Output File | Description |
|------|------------|-------------|
| Step 1 | `variants.vcf` | Raw variant calls |
| Step 2 | `filtered_variants.vcf` | Quality-filtered variants |
| Step 3 | `filtered_variants_annotated.vcf` | Functional impact analysis |
| Step 4 | `high_impact_variants.vcf` | Severe impact variants |
| Step 5 | `filtered_variants_annotated_with_dbSNP.vcf` | Matched variants with dbSNP |
| Step 6 | `*.sorted.bam`, `*.bai` | Indexed BAM files |

---
### Summary of Variant Analysis

#### Variant Types:
SNPs: 191,547 (99.92%)
Indels: 149 (0.08%)

#### Quality Scores:
Median: 199.503
Range: 3.03 to 3549.08

#### Mutation Patterns:
Transitions: 182,589 (95.14%)
Transversions: 9,319 (4.86%)
Ti/Tv Ratio: 19.59


