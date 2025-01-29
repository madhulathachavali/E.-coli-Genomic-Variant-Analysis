#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam = "results/alignment/aligned/*.sorted.bam"
params.ref_genome = "/home/akella/madhu/Ecoli/ASM584v2/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
params.outdir = "results/variant_calling"

process VARIANT_CALLING {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_genome

    output:
    path "${sample_id}.vcf"

    script:
    """
    bcftools mpileup -Ou -f ${reference_genome} ${bam} | \
    bcftools call --ploidy 1 -mv -Ov -o ${sample_id}.vcf
    """
}

workflow {
    bam_ch = channel
        .fromPath(params.bam)
        .map { file -> 
            def sample_id = file.name.replaceAll(/\.sorted\.bam$/, '')
            tuple(sample_id, file, file + ".bai")
        }
    reference_ch = channel.fromPath(params.ref_genome)

    VARIANT_CALLING(bam_ch, reference_ch)
}

