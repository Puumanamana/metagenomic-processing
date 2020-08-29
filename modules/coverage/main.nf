nextflow.enable.dsl = 2

process bwa_index {
    container 'nakor/metagenome-assembly'
    
    input:
    tuple val(label), path(assembly)

    output:
    tuple val(label), path('ref*')

    script:
    """
    bwa index ${assembly} -p ref
    """
}

process bwa {
    publishDir params.outdir, mode: 'copy'
    container 'nakor/metagenome-assembly'
    
    input:
    tuple val(sample), path(fastq), path(index)

    output:
    tuple val(sample), path('coverage*.bam'), emit: bam
    tuple val(sample), path('coverage*.bai'), emit: bai
    
    script:
    """
    bwa mem -a -M -t $task.cpus ref $fastq \
        | samtools sort -@ $task.cpus -o coverage_${sample}.bam
    samtools index -@ $task.cpus coverage_${sample}.bam coverage_${sample}.bai
    """    
}

process aln_stats {
    publishDir params.outdir, mode: 'copy'
    container 'nakor/metagenome-assembly'
    
    input:    
    tuple val(sample), path(bam), path(bai), path(fasta)
    each minLen
    each flag
    each qual

    output:
    stdout

    script:
    """
    #!/usr/bin/env bash
    
    # -f 3 includes reads mapped in proper pair
    # -F 2308 excludes supplementary alignments and unmapped reads

    bioawk -c fastx '{OFS="\\t"}{if (length(\$seq)>$minLen) print \$name,1,length(\$seq)}' $fasta > ctg_list.bed
    count=\$(samtools view -c -f ${flag[1]} -F 2308 -q $qual -@ $task.cpus $bam -L ctg_list.bed)
    
    echo "bwa,>$minLen bp,$qual,${flag[1]},$sample,\$count" | tr -d '\\n'
    """
}

workflow coverage {
    take:
    reads
    assembly // assembly needs to be prefixed with sample name if not co-assembly

    main:
    assembly | unique{ it[0] } | bwa_index

    if(params.coassembly) {
        reads.combine(bwa_index.out).map{[it[0], it[1], it[3]]} | bwa
        cov_and_assembly = bwa.out.bam.join(bwa.out.bai).combine(assembly.map{it[1]})
    } else {
        reads.combine(bwa_index.out, by: 0) | bwa
        cov_and_assembly = bwa.out.bam.join(bwa.out.bai).combine(assembly, by: 0)
    }

    aln_stats(cov_and_assembly, params.cov_summary_len, params.cov_summary_flags,
              params.cov_summary_quals)
        
    emit:
    all = cov_and_assembly
    bam = bwa.out.bam
    bai = bwa.out.bai
    summary = aln_stats.out
}

workflow {
    reads = Channel.fromFilePairs(params.reads)
    assembly = Channel.fromPath(params.assembly)
    coverage(reads, assembly)
}

workflow test {
    reads = Channel.fromFilePairs("$baseDir/../../test_data/*_R{1,2}.fastq.gz")
    assembly = Channel.fromPath("$baseDir/../../test_data/assembly_sample*.fasta")
        .map{[it.getSimpleName().tokenize('_')[1], it]}
    coverage(reads, assembly)
}
