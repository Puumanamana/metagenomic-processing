nextflow.enable.dsl=2

// Default cropping values
// params.hd_crop ?: 0
// params.tl_crop ?: 0
// params.hd_crop_fwd ?: params.hd_crop
// params.hd_crop_rev ?: params.hd_crop
// params.tl_crop_fwd ?: params.tl_crop
// params.tl_crop_rev ?: params.tl_crop

include { reads_qc; trimming } from './modules/trimming.nf' addParams(outdir: "${params.outdir}/preprocessing")
include { assembly } from './modules/assembly.nf' addParams(outdir: "${params.outdir}/assembly")
include { coverage; aln_stats } from './modules/coverage.nf' addParams(outdir: "${params.outdir}/coverage")
include { annotation } from './modules/annotation.nf' addParams(outdir: "${params.outdir}/annotation")
include { binning } from './modules/binning.nf' addParams(outdir: "${params.outdir}/binning")

process reshape_summary {
    publishDir params.outdir

    input:
    file(summary)

    output:
    file('summary_table.csv')

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    table = pd.read_csv("${summary}", index_col=0, header=None, names=['step', 'sample', 'count'])
    table = table.pivot('sample', 'step', 'count').reindex(columns=table.step.unique())

    table.to_csv('summary_table.csv')
    """
}

workflow wgs_analysis {
    take: reads
    main:
    reads | reads_qc
    reads | trimming | assembly

    if(params.split_assembly) {
        assembly.out | join(trimming.out) | map{[it[1], it[0], it[2]]} | coverage
        coverage.out | join(assembly.out) | set{coverage_and_assembly}
    } else {
        assembly.out | map{it[1]} | combine(trimming.out) | coverage
        coverage.out | combine(assembly.out.map{it[1]}) | set{coverage_and_assembly}
        assembly.out | map{it[1]} | combine(coverage.out.map{it[1..-1]}) | groupTuple | binning
    }
    assembly.out | map{it[1]} | annotation

    // Counts summary
    reads | map{"0,raw,${it[0]},${it[1][0].countFastq()}"} | set{raw_seq_counts}
    trimming.out | map{"1,trimming,${it[0]},${it[1][0].countFastq()}"} | set{trim_seq_counts}
    assembly.out | map{"2,${params.assembler},${it[0]},${it[1].countFasta()}"} | set{assembly_counts}
    coverage_and_assembly | combine(Channel.from(0, 500, 1000, 2000, 5000)) | aln_stats

    
    raw_seq_counts| mix(trim_seq_counts, assembly_counts, aln_stats.out) \
        | collectFile(newLine: true, sort: true) \
        | reshape_summary

}

workflow {    
    Channel.fromFilePairs(params.reads) | wgs_analysis
}
