#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "$projectDir/results/trimmed/*_{1,2}.fastq"
params.ref_genome = "/home/akella/madhu/Ecoli/ASM584v2/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
params.outdir = "results"

process INDEX_GENOME {
    publishDir "${params.outdir}/bwa_index", mode: 'copy'

    input:
    path ref_genome

    output:
    path "${ref_genome}*"

    script:
    """
    bwa index ${ref_genome}
    """
}

process ALIGN_READS {
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem ${index[0]} ${reads} | samtools view -bS - > ${sample_id}.bam
    """
}

workflow {
    reference_ch = channel.fromPath(params.ref_genome)
    reads_ch = channel.fromFilePairs(params.reads)

    index = INDEX_GENOME(reference_ch)
    ALIGN_READS(reads_ch, index)
}

