spades_modes = [['spades', ''], ['metaspades', '--meta']]
input_modes = ['paired-only', 'all']

Channel.fromFilePairs(params.reads)
    .multiMap{it ->
	qc: it
	trimming: it}
    .set{RAW_FASTQ}

process PreQC {
    tag {"preQC-${sample}"}
    publishDir params.outdir+"/0-preQC", mode: "copy"

    input:
    tuple val(sample), file(fastqs) from RAW_FASTQ.qc

    output:
    file("*.html")
    file("*.zip") into PRE_QC_FOR_MULTIQC

    script:
    """
    fastqc -t ${task.cpus} *.fastq*
    """
}

process PreMultiQC {
    tag {"preMultiQC"}
    publishDir params.outdir+"/0-preQC", mode: "copy"

    input:
    file f from PRE_QC_FOR_MULTIQC.collect()

    output:
    file("multiqc*.html")

    script:
    """
    multiqc .
    """
}

process Trimming {
    // Quality filter and trimming
    tag { "trimmomatic-${sample}" }
    publishDir params.outdir+"/1-trimming", mode: "copy"

    input:
    tuple val(sample), file(fastqs) from RAW_FASTQ.trimming

    output:
    tuple val(sample), file("*paired_R*.fastq.gz") into TRIMMED_FASTQ

    script:
    """
    #!/usr/bin/env bash

    # [ "${params.adapters}" = "null" ] && args="" || args="ILLUMINACLIP:${params.adapters}:2:30:10:2:keepBothReads"
    args=""

    java -jar ${HOME}/.local/src/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} \
        ${fastqs} \
        ${sample}_paired_R1.fastq.gz ${sample}_unpaired_R1.fastq.gz \
	${sample}_paired_R2.fastq.gz ${sample}_unpaired_R2.fastq.gz \
	\$args LEADING:3 MINLEN:100 
    """
}

TRIMMED_FASTQ.multiMap{it ->
    qc: it
    assembly: it[1]
    coverage: it
}.set{FILTERED_FASTQ}

process PostQC {
    tag {"postQC_${sample}"}
    publishDir params.outdir+"/2-postQC", mode: "copy"

    input:
    tuple val(sample), file(fastqs) from FILTERED_FASTQ.qc

    output:
    file("*.html")
    file("*.zip") into POST_QC_FOR_MULTIQC

    script:
    """
    fastqc *fastq.gz
    """
}


process PostMultiQC {
    tag {"postMultiQC"}
    publishDir params.outdir+"/2-postQC", mode: "copy"

    input:
    file f from POST_QC_FOR_MULTIQC.collect()

    output:
    file("multiqc*.html")

    script:
    """
    multiqc .
    """
}

process MetaspadesAssemby {
    tag {"metaSPAdes"}
    publishDir params.outdir+"/3-assemblies", mode: "copy"

    input:
    each assembly_mode from spades_modes
    each input_mode from input_modes
    file f from FILTERED_FASTQ.assembly.collect()

    output:
    tuple(val(assembly_mode), val(input_mode), file("assembly_*.fasta")) into ASSEMBLY
    
    script:
    """
    cat *_paired_R1.fastq.gz > all_R1.fastq.gz
    cat *_paired_R2.fastq.gz > all_R2.fastq.gz
    cat *_unpaired*.fastq.gz > all_unpaired.fastq.gz

    spades_path="${HOME}/.local/src/SPAdes-3.14.1-Linux/bin"
    
    [ "${input_mode}" = "paired-only" ] && input_mode="" || input_mode="-s all_unpaired.fastq.gz"
    \${spades_path}/spades.py ${assembly_mode[1]} -k 21,33,55,77 \
        -1 all_R1.fastq.gz -2 all_R2.fastq.gz \$input_mode \
        -t ${task.cpus} -m 80 -o spades_output

    mv spades_output/scaffolds.fasta "assembly_${assembly_mode[0]}_${input_mode}.fasta"
    """    
}

process CoverageIndex {
    tag {"coverageIndex-${sample}-${assembly_mode[0]}-${input_mode}"}

    input:
    tuple(val(assembly_mode), val(input_mode), file(assembly)) from ASSEMBLY

    output:
    tuple(val(assembly_mode), val(input_mode), file('db*')) into BWA_DB

    script:
    """
    bwa index ${assembly} -p db
    """
	  
}

process Coverage {
    tag {"coverage_${sample}-${assembly_mode[0]}-${input_mode}"}
    publishDir params.outdir+"/4-coverage", mode: "copy"

    input:
    tuple(val(sample), file(trimmed_reads), val(assembly_mode), val(input_mode), file(index)) from FILTERED_FASTQ.coverage.combine(BWA_DB)

    output:
    file('coverage*.bam') into COVERAGE 
    
    script:
    """
    bwa mem -a -M -t 50 db ${trimmed_reads} \
        | samtools sort -@ ${task.cpus} -o coverage_${assembly_mode[0]}_${input_mode}_${sample}.bam
    """    
}

process Binning {
}
