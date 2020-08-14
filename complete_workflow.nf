nextflow.enable.dsl=2

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
    table = table.applymap(lambda x: "{:,}k".format(int(x/1000)))

    table.to_csv('summary_table.csv')
    """
}

workflow wgs_analysis {
    take: reads
    main:
    reads | reads_qc
    reads | trimming | assembly

    if(params.split_assembly) {
        assembly.out | combine(trimming.out, by: 0) | map{[it[1], it[0], it[2]]} | coverage
        coverage.out | combine(assembly.out, by: 0) | set{coverage_and_assembly}
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

    aln_stats(coverage_and_assembly,
              [0, 1000, 2000, 5000], // min ctg lengths
              ['mapped', 0], ['mapped in proper pair', 3], // flags
              [0, 30, 50] // min mapping qualities
    )
    
    raw_seq_counts| mix(trim_seq_counts, assembly_counts, aln_stats.out) \
        | collectFile(newLine: true, sort: true) \
        | reshape_summary

}

workflow {
    Channel.fromFilePairs(params.reads) | wgs_analysis
}

workflow test {
    Channel.fromFilePairs("$baseDir/test_data/sample*_R{1,2}.fastq.gz") | wgs_analysis
}
