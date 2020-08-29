nextflow.enable.dsl = 2

process coconet {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/coconet", mode: 'copy'
    container 'nakor/metagenome-assembly'
    
    input:
    file(fasta)
    file(bam)
    file(bai)

    output:
    file('coconet_bins*.csv') 

    script:
    """
    coconet --fasta $fasta --coverage $bam --output output
    mv output/bins_*.csv coconet_bins.csv
    """
}

process metabat2 {
    publishDir "${params.outdir}/metabat2", mode: 'copy'
    container 'nakor/metagenome-assembly'
    
    input:
    file(fasta)
    file(bam)
    file(bai)
    
    output:
    file('metabat2_bins*.csv')

    script:
    """
    runMetaBat.sh --saveCls $fasta $bam
    mv *metabat-bins*/bin metabat2_bins.csv
    """
}

workflow binning {
    take:
    fasta
    bam
    bai
    
    main:
    metabat2(fasta, bam, bai)
    coconet(fasta, bam, bai)    
}

workflow {
    fasta = Channel.fromPath(params.assembly)
    bam = Channel.fromPath(params.coverage)
    bai = Channel.fromPath(params.coverage_index)

    binning(fasta, bam.collect(), bai.collect())
}

workflow test {
    fasta = Channel.fromPath("$baseDir/../../test_data/assembly.fasta")
    bam = Channel.fromPath("$baseDir/../../test_data/*.bam")
    bai = Channel.fromPath("$baseDir/../../test_data/*.bai")

    binning(fasta, bam.collect(), bai.collect())
}

