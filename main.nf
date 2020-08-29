nextflow.enable.dsl=2

include { reads_qc; trimming } from './modules/trimming/main.nf' addParams(outdir: "${params.outdir}/preprocessing")
include { assembly } from './modules/assembly/main.nf' addParams(outdir: "${params.outdir}/assembly")
include { coverage; aln_stats } from './modules/coverage/main.nf' addParams(outdir: "${params.outdir}/coverage")
include { annotation } from './modules/annotation/main.nf' addParams(outdir: "${params.outdir}/annotation")
include { binning } from './modules/binning/main.nf' addParams(outdir: "${params.outdir}/binning")
include { reshape_summary } from './modules/misc/main.nf'

workflow wgs_analysis {
    take: reads
    main:
    reads | reads_qc
    reads | trimming | assembly

    coverage(trimming.out, assembly.out)

    if (params.coassembly) {
        binning(assembly.out, coverage.out.bam.collect(), coverage.out.bai.collect())
    }

    assembly.out | map{it[1]} | annotation

    // Counts summary
    reads | map{"0,raw,${it[0]},${it[1][0].countFastq()}"} | set{raw_seq_counts}
    trimming.out | map{"1,trimming,${it[0]},${it[1][0].countFastq()}"} | set{trim_seq_counts}
    assembly.out | map{"2,${params.assembler},${it[0]},${it[1].countFasta()}"} | set{assembly_counts}

    raw_seq_counts| mix(trim_seq_counts, assembly_counts) \
        | collectFile(newLine: true, sort: true, storeDir: params.outdir,
                      name: 'preproc_summary.csv')
    coverage.out.summary.collectFile(newLine: true, sort: true, storeDir: params.outdir,
                                 name: 'coverage_summary.csv')

}

workflow {
    Channel.fromFilePairs(params.reads) | wgs_analysis
}

workflow test {
    Channel.fromFilePairs("$baseDir/test_data/sample*_R{1,2}.fastq.gz") | wgs_analysis
}
