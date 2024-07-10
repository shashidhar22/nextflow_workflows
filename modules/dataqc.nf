#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Run FASTQC
//
// PROCESS INPUTS:
// - Tuple with sample name and path to FASTQ files
// - Output path based on the stage of analysis (pre-processing/trimmed)
// PROCESS OUTPUTS:
// - FASTQC zip and html outputs
// EMITTED CHANNELS:
// - fastqc_results :- FASTQC ZIP reports
// - fastqc_reports :- FASTQC HTML reports
//
// NOTE: 
// It is recommended not to chain QC modules, a review will be required at the 
// end of pre-processing steps. 
// TODO: 
process runFastQC {
  conda "$params.qc_env"
  label 'low_mem'
  publishDir "$params.output.data/preprocess/fastqc/${mode}/archives/", 
    mode : "copy", pattern : "*.zip"
  publishDir "$params.output.data/preprocess/fastqc/${mode}/html",
    mode : "copy", pattern : "*.html"
  input:
    tuple val(sample), path(fastq_paths)
    val mode
  output:
    path "*.zip", emit: fastqc_results
    path "*.html", emit: fastqc_reports
  script:
    """
    fastqc -t $params.qc_threads --memory $params.qc_mem \
      --noextract -q  -o ./ *.f*q*  > stdout.txt 2> stderr.txt
    """
}


// Run MULTIQC
//
// PROCESS INPUTS:
// - Outputs from each of the QC stages
// PROCESS OUTPUTS:
// - MULTIQC reports
// EMITTED CHANNELS:
//
// NOTE: 
// It is recommended not to chain QC modules, a review will be required at the 
// end of pre-processing steps. 
// TODO: 
process runMultiQC {
  cache false
  conda "$params.qc_env"
  label 'mid_mem'
  publishDir "$params.output.data/preprocess/multiqc/", mode : "copy"

  input:
    path pre_fastqc_results
    path bbduk_results
    path post_fastqc_results

  output:
    path "*", emit: multiqc_results
    
  script:
    """
    python -m multiqc -c "$params.qc_report_config" -i "$params.output.qc_report_title"   .
    """
}



// Trim reads using BBDuk to remove adapter sequences and lower quality bases
//
// PROCESS INPUTS:
// - Tuple with sample name, path to raw FASTQ file, and path to adapter file
// (Source: fastq_path.concat(generateFastq.out.sim_fastq).combine(adapter_path))
//
// PROCESS OUTPUTS:
// - Trimmed FASTQ files with adapter and low quality bases removed
// EMITTED CHANNELS:
// - trimmed_reads :- Channel containing sample name, and path to trimmed FASTQ 
// file
//
// NOTE: 
//
// TODO: 
process runBBduk {
  conda "$params.qc_env"
  label 'mid_mem'
  publishDir "$params.output.data/preprocess/trimmed_reads", mode : "copy"
  input:
    tuple val(sample), path(fastq_path)
  output:
    path "*clean.fq", emit: trimmed_reads
    tuple val(sample), path("*clean.fq"), emit: trimmed_qc
    path "${sample}_bbduk_metrics.log", emit: trimmed_stats
  script:
  if ( fastq_path.size() == 2 )
    """
    bbduk.sh -Xmx4g in1=${fastq_path[0]} in2=${fastq_path[1]} \
    out1=${sample}_R1.clean.fq out2=${sample}_R2.clean.fq \
    ref=adapters qtrim=rl stats=${sample}_bbduk_stats.txt \
    ktrim=r k=23 mink=11 hdist=1 trimq=30 > stdout.log \
    2> ${sample}_bbduk_metrics.log
    """
  else 
    """
    bbduk.sh -Xmx4g in1=${fastq_path[0]} \
    out1=${sample}.clean.fq stats=${sample}_bbduk_stats.txt\
    ref=adapters qtrim=rl ktrim=r k=23 mink=11 hdist=1 trimq=30 \
    > stdout.log 2> ${sample}_bbduk_metrics.log
    """
} 
