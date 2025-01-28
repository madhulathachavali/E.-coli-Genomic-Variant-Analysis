# E.-coli-Genomic-Variant-Analysis

### 
This project focuses on processing 15 E. coli samples to generate and analyze VCF files using Nextflow. The workflow uses command-line tools such as Nextflow, bcftools, and the SRA Toolkit

### Objectives
- Perform quality control (QC) on raw sequencing data.
- Convert raw data to SAM, BAM, and eventually VCF format.
- Perform variant calling, allele frequency analysis, and phylogenetic studies.

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
