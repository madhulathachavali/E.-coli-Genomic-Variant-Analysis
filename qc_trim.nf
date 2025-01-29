#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "$projectDir/fastq_output/*_{1,2}.fastq"
params.single_end_reads = "$projectDir/fastq_output/*_1.fastq"
params.outdir = "results"
params.trimmomatic = "/usr/share/java/trimmomatic-0.39.jar"

process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process TRIMMOMATIC {
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed*.fastq")

    script:
    if (reads instanceof List && reads.size() == 2)
        """
        java -jar ${params.trimmomatic} PE -threads 8 -phred33 \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_trimmed_1.fastq ${sample_id}_unpaired_1.fastq \
        ${sample_id}_trimmed_2.fastq ${sample_id}_unpaired_2.fastq \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    else
        """
        java -jar ${params.trimmomatic} SE -threads 8 -phred33 \
        ${reads} \
        ${sample_id}_trimmed.fastq \
        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
}

workflow {
    Channel
        .fromFilePairs(params.reads, size: -1)
        .mix(Channel.fromPath(params.single_end_reads).map { file -> tuple(file.simpleName, file) })
        .set { read_pairs_ch }

    FASTQC(read_pairs_ch)
    TRIMMOMATIC(read_pairs_ch)
}



