nextflow.enable.dsl=2

params.reads = "/home/akella/madhu/Ecoli/fastq_output/*_{1,2}.fastq"
params.trimmed_reads = "/home/akella/madhu/Ecoli/results/trimmed/*_trimmed_{1,2}.fastq"
params.outdir = "results/multiqc"

process FASTQC_RAW {
    publishDir "${params.outdir}/raw_fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.zip"
    
    script:
    """
    fastqc ${reads}
    """
}

process FASTQC_TRIMMED {
    publishDir "${params.outdir}/trimmed_fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.zip"
    
    script:
    """
    fastqc ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path raw_fastqc
    path trimmed_fastqc
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc ${raw_fastqc} ${trimmed_fastqc} -o .
    """
}

workflow {
    raw_reads_ch = channel.fromFilePairs(params.reads)
    trimmed_reads_ch = channel.fromFilePairs(params.trimmed_reads)
    
    raw_fastqc = FASTQC_RAW(raw_reads_ch)
    trimmed_fastqc = FASTQC_TRIMMED(trimmed_reads_ch)
    
    MULTIQC(raw_fastqc.collect(), trimmed_fastqc.collect())
}

