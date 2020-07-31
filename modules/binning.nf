nextflow.enable.dsl = 2

process coconet {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/coconet", mode: 'copy'

    input:
    tuple(file(fasta), file(bam), file(bai))

    output:
    file('coconet_bins*.csv') 

    script:
    """
    coconet --fasta $fasta --coverage $bam --output output
    mv output/bins_*.csv coconet_bins.csv
    """
}

process metabat2 {
    container "metabat/metabat"
    publishDir "${params.outdir}/metabat2", mode: 'copy'

    input:
    tuple(file(fasta), file(bam), file(bai))

    output:
    file('metabat2_bins*.csv')

    script:
    """
    runMetaBat.sh --saveCls $fasta $bam
    mv *metabat-bins*/bin metabat2_bins.csv
    """
}

workflow binning {
    take: data
    main: data | (metabat2 & coconet)
}

workflow {
    fasta = Channel.fromPath(params.assembly)
    coverage = Channel.fromPath(params.coverage).collect()

    fasta | combine(coverage) | map{[it[0], it[1..-1]]} | binning
}
