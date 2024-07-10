#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Build Index for reference genomes 
//
// PROCESS INPUTS:
// - Reference FASTA file path
// - GTF file
// PROCESS OUTPUTS:
// - Genome index
// EMITTED CHANNELS:
// - star_index :- Genome index
//
// NOTE: 
// This module will only need to be run in the absence of an index
// TODO: 
process buildSTARIndex {
  conda "$params.align_env"
  label 'high_mem'
  publishDir "$params.output.reference/index/star/${mode}/${genome}/", 
    mode : "copy"
  input:
    each reference_param
    each star_param
  output:
    path "*", emit: star_index
  script:
    genome = reference_param[0]
    reference_fasta = reference_param[1]
    reference_gtf = reference_param[2]
    star_sanbases = reference_param[3]
    mode = star_param[0]
    star_overhang = star_param[1]

    """
    STAR --runThreadN $params.index_threads --runMode genomeGenerate \
     --genomeDir ./ --genomeFastaFiles ${reference_fasta} \
     --genomeSAindexNbases ${star_sanbases} \
     --sjdbGTFfile ${reference_gtf} \
     --sjdbOverhang ${star_overhang}
    """
}
