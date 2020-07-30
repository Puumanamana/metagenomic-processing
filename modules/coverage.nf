nextflow.enable.dsl = 2

process bwa_index {
    input:
    tuple(val(label), file(assembly))

    output:
    tuple(val(label), file('ref*'))

    script:
    """
    bwa index ${assembly} -p ref
    """
}

process bwa {
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple(val(sample), file(fq), file(index))

    output:
    file('coverage*.bam')
    
    script:
    """
    bwa mem -a -M -t 50 ref $fq \
        | samtools sort -@ $task.cpus -o coverage_${sample}.bam
    """    
}

process bowtie2_index {
    input:
    tuple(val(label), file(assembly))

    output:
    tuple(val(label), file('ref*'))

    script:
    """
    bowtie2-build --threads $task.cpus $assembly ref
    """
}

process bowtie2 {
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple(val(sample), file(fq), file(index))

    output:
    file('coverage*.bam')
    
    script:
    """
    bowtie2 -p $task.cpus -x ref -1 ${fq[0]} -2 ${fq[1]} | \
            samtools view -@ $task.cpus -bS | \
            samtools sort -@ $task.cpus -o ${sample}.bam
    samtools index -@ "$task.cpus ${sample}.bam ${sample}.bai
    """    
}

workflow coverage {
    // data is a (assembly, sample, fastq tuple) channel
    take: data

    main:
    data | multiMap {it ->
        assembly: [it[0].getSimpleName(), it[0]]
        reads: [it[0].getSimpleName(), it[1], it[2]]
    } | set{split}
    
    if (params.aligner == "bwa") {
        split.assembly | unique{ it[0] } | bwa_index | set{reference}
        split.reads | combine(reference, by: 0) | map{it[1..-1]} | bwa | set{bamcov}
    }
    else {
        split.assembly | unique{ it[0] } | bowtie2_index | set{reference}
        split.reads | combine(reference, by: 0) | map{it[1..-1]} | bowtie2 | set{bamcov}
    }
    emit: bamcov
}

workflow {
    reads = Channel.fromFilePairs(params.reads)
    assembly = Channel.fromPath(params.assembly)
    assembly | combine(reads) | coverage
}
