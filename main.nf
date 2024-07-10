nextflow.enable.dsl=2
include {runFastQC as preFastQC ; runFastQC as postFastQC ; runMultiQC ; runBBduk} from './modules/dataqc'
include { buildSTARIndex } from './modules/alignment'
// FASTQ data QC
workflow data_qc {
  fastq_path = Channel.fromFilePairs(params.data.fastq_path)
  main:
    // Run FASTQC on all input FASTQ files
    preFastQC(fastq_path, "raw")
    //Run BBDuk to trim all low quality and adapter reads
    runBBduk(fastq_path)
    // RunFASTQC again trimmed reads
    postFastQC(runBBduk.out.trimmed_qc, "trimmed")
    //Run MULTIQC report
    runMultiQC(preFastQC.out.fastqc_results.collect(),
        runBBduk.out.trimmed_stats.collect(),
        postFastQC.out.fastqc_results.collect())
    
}

// Build Alignment indices for Human, HIV, KSHV
workflow build_index {
  organism = Channel.from(params.reference.organism)
  reference_fasta_path = Channel.fromPath(params.reference.fasta)
  reference_gtf_path = Channel.fromPath(params.reference.gtf)
  reference_fasta = organism.merge(reference_fasta_path)
  reference_gtf = organism.merge(reference_gtf_path)
  star_sanbases = Channel.from(params.star.index.sanbases)
  star_param = Channel.from(params.star.index.overhang)
  reference_params = reference_fasta.join(reference_gtf.join(star_sanbases))
  main:
    // Run STAR index
    buildSTARIndex(reference_params, star_param)
}