#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define sample IDs explicitly
params.samples = [
    'SRR13921543', 'SRR13921545', 'SRR13921546', 
    'SRR13921549', 'SRR13921550', 'SRR13921551', 
    'SRR13921552', 'SRR13921554', 'SRR13921555', 
    'SRR13921556', 'SRR13921557', 'SRR13921558', 
    'SRR13921559', 'SRR13921560', 'SRR13921561'
]

params.reads_dir = "$HOME/madhu/Ecoli/results/processed_samples/fastq/paired_end"
params.ref_genome = "/home/akella/madhu/Ecoli/ASM584v2/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
params.outdir = "results/aligned"

workflow {
    // Create channel with sample IDs and their corresponding read files
    read_pairs_ch = Channel
        .fromList(params.samples)
        .map { sample_id -> 
            def read1 = "${params.reads_dir}/${sample_id}_trimmed_1.fastq"
            def read2 = "${params.reads_dir}/${sample_id}_trimmed_2.fastq"
            [sample_id, file(read1), file(read2)]
        }
        .view { sample -> "Processing Sample: ${sample[0]}" }
}

