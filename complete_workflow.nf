nextflow.enable.dsl=2

include { reads_qc; trimming } from './modules/trimming.nf' addParams(outdir: "${params.outdir}/preprocessing")
include { assembly } from './modules/assembly.nf' addParams(outdir: "${params.outdir}/assembly")
include { coverage } from './modules/coverage.nf' addParams(outdir: "${params.outdir}/coverage")
include { annotation } from './modules/annotation.nf' addParams(outdir: "${params.outdir}/annotation")
include { binning } from './modules/binning.nf' addParams(outdir: "${params.outdir}/binning")

workflow wgs_analysis {
    take: reads
    main:
    reads | reads_qc
    reads | trimming | assembly

    if(params.split_assembly) {
        assembly.out | combine(trimming.out, by: 0) | map{it[1..-1]} | coverage
    } else {
        assembly.out | map{it[1]} | combine(trimming.out) | coverage
        assembly.out | map{it[1]} | combine(coverage.out) | groupTuple | binning
    }
    assembly.out | map{it[1]} | annotation
}

workflow {
    Channel.fromFilePairs(params.reads) | wgs_analysis
}
