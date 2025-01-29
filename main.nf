params.sra_dir = "/home/akella/madhu/Ecoli/sra"
params.outdir = "fastq_output"

process FASTQ_DUMP {
    publishDir params.outdir, mode: 'copy'

    input:
    path sra_file

    output:
    path "${sra_file.baseName}*.fastq"

    script:
    """
    fastq-dump ${sra_file} -O . --split-files
    """
}
workflow {
    Channel
        .fromPath("${params.sra_dir}/*.sra")
        .set { sra_files }

    FASTQ_DUMP(sra_files)
}

