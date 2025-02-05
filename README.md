# E.-coli-Genomic-Variant-Analysis

### 
This project focuses on processing 15 E. coli samples to generate and analyze VCF files using Nextflow. The workflow uses command-line tools such as Nextflow, bcftools, and the SRA Toolkit (https://www.ncbi.nlm.nih.gov/sra?term=SAMN02604091)

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
- bcftools: Variant calling and VCF file generation.
- Phylogenetic Analysis Tools: For studying evolutionary relationships.

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



