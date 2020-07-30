nextflow.enable.dsl = 2

process trimmomatic {
    // Quality filter and trimming
    publishDir "${params.outdir}/trimming", mode: 'copy'
    
    input:
    tuple(val(sample), file(fastqs))
    
    output:
    tuple val(sample), file("*.fastq.gz")

    script:
    """
    #!/usr/bin/env bash

    [ "${params.adapters}" = "null" ] && opt_args="" \
        || opt_args="ILLUMINACLIP:${params.adapters}:2:30:10:2:keepBothReads"

    avg_read_len=\$(zcat ${fastqs[0]} | head -n 4000 | awk '{if(NR%4==2) {count++; bases += length} } END{print int(bases/count)}')
    tl_crop=\$((\$avg_read_len-${params.tl_crop}))

    java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} \
         ${fastqs} \
         ${sample}_paired_R1.fastq.gz ${sample}_unpaired_R1.fastq.gz \
         ${sample}_paired_R2.fastq.gz ${sample}_unpaired_R2.fastq.gz \
         \$opt_args MINLEN:${params.min_trim_len} \
         AVGQUAL:${params.min_avg_qual} \
         LEADING:${params.clip_qual} TRAILING:${params.clip_qual} \
         SLIDINGWINDOW:${params.trim_win_len}:${params.trim_win_qual} \
         HEADCROP:${params.hd_crop} CROP:\$tl_crop

    [ ${params.keep_unpaired} = false ] && rm *unpaired* || echo "Keeping unpaired reads"\
    """
}

process fastp {
    // Quality filter and trimming
    publishDir "${params.outdir}/trimming", mode: 'copy'
    
    input:
    tuple(val(sample), file(fastqs))

    output:
    tuple val(sample), file("*.fastq.gz")
    
    script:
    """
    #!/usr/bin/env bash

    [ "${params.adapters}" = "null" ] && adapter_arg="" \
        || adapter_arg="--adapter_fasta ${params.adapters}"

    fastp --trim_poly_g -w ${task.cpus} --average_qual ${params.min_avg_qual} \
          --length_required ${params.min_trim_len} --cut_front --cut_tail \
          --cut_mean_quality ${params.trim_win_qual} --cut_window_size ${params.trim_win_len} \
          --trim_front1 ${params.hd_crop} --trim_tail1 ${params.tl_crop} \$adapter_arg \
          -i ${fastqs[0]} -I ${fastqs[1]} \
          -o ${sample}_paired_R1.fastq.gz -O ${sample}_paired_R2.fastq.gz \
          --unpaired1 ${sample}_unpaired_R1.fastq.gz \
          --unpaired2 ${sample}_unpaired_R2.fastq.gz

    mv fastp.html fastp_${sample}.html
    """
}

process fastqc {
    publishDir "${params.outdir}/readsQC/fastQC", pattern: "*.html", mode: 'copy'
    
    input:
    tuple val(sample), file(fastq)
    
    output:
    file("*.{zip,html}")

    script:
    """
    fastqc -t ${task.cpus} *.fastq*
    """
}

process multiqc {
    publishDir "${params.outdir}/readsQC/multiQC", pattern: "*.html", mode: 'move'
    
    input:
    file(fastqs)
    
    output:
    file("multiqc.html")

    script:
    """
    multiqc . -n multiqc.html
    """
}

// Workflows
workflow reads_qc {
    take: data
    main: data | fastqc | collect | multiqc
}

workflow trimming {
    // data is a (sample, [fwd, rev]) channel
    take: data
    
    main:
    if (params.trimming == "fastp") {data | fastp | set{trimmed_reads}}
    else {data | trimmomatic | set{trimmed_reads}}

    if (!params.keep_unpaired) {
        trimmed_reads | transpose | filter{ !(it[1].baseName =~/.*unpaired.*/) } \
            | groupTuple \
            | set{filtered_reads}
    } else {
        trimmed_reads | set{filtered_reads}
    }
    
    filtered_reads | reads_qc
    
    emit: filtered_reads
}

workflow {
    fastqs = Channel.fromFilePairs(params.reads)
    fastqs | (reads_qc & trimming)
}
