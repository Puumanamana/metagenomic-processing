spades_modes = [['metaspades', '--meta']]
input_modes = ['paired-only']

Channel.fromFilePairs(params.reads)
    .multiMap{it ->
    qc: it
    trimming: it}
    .set{RAW_FASTQ}

process PreQC {
    tag {"preQC-${sample}"}
    publishDir params.outdir+"/0-preQC", mode: "copy", pattern: "*.html"

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

if (params.trimming == 'trimmomatic') {
    process Trimmomatic {
        // Quality filter and trimming
        tag { "trimmomatic-${sample}" }
        publishDir params.outdir+"/1-trimming", mode: "copy"

        input:
        tuple val(sample), file(fastqs) from RAW_FASTQ.trimming

        output:
        tuple val(sample), file("*_paired_R*.fastq.gz"), file("*_unpaired.fastq.gz") into TRIMMED_FASTQ

        script:
        """
        #!/usr/bin/env bash

        [ "${params.adapters}" = "null" ] && args="" || args="ILLUMINACLIP:${params.adapters}:2:30:10:2:keepBothReads"
        args=""

        avg_read_len=\$(zcat ${fastqs[0]} | head -n 4000 | awk '{if(NR%4==2) {count++; bases += length} } END{print int(bases/count)}')
        tl_crop=\$((\$avg_read_len-${params.tl_crop}))

        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} \
             ${fastqs} \
             ${sample}_paired_R1.fastq.gz ${sample}_unpaired_R1.fastq.gz \
             ${sample}_paired_R2.fastq.gz ${sample}_unpaired_R2.fastq.gz \
             \$args MINLEN:${params.min_trim_len} \
             AVGQUAL:${params.min_avg_qual} \
             LEADING:${params.clip_qual} TRAILING:${params.clip_qual} \
             SLIDINGWINDOW:${params.trim_win_len}:${params.trim_win_qual} \
             HEADCROP:${params.hd_crop} CROP:\$tl_crop

        cat *_unpaired_R*.fastq.gz > ${sample}_unpaired.fastq.gz
        """
    }
} else if (params.trimming == 'fastp'){
    process Fastp {
        // Quality filter and trimming
        tag { "fastp-${sample}" }
        publishDir params.outdir+"/1-trimming", mode: "copy"

        input:
        tuple val(sample), file(fastqs) from RAW_FASTQ.trimming

        output:
        tuple val(sample), file("*_paired_R*.fastq.gz"), file("*_unpaired.fastq.gz") into TRIMMED_FASTQ
        file("*.html")

        script:
        """
        #!/usr/bin/env bash

        fastp --trim_poly_g -w ${task.cpus} --average_qual ${params.min_avg_qual} \
              --length_required ${params.min_trim_len} --cut_front --cut_tail \
              --cut_mean_quality ${params.trim_win_qual} --cut_window_size ${params.trim_win_len} \
              --trim_front1 ${params.hd_crop} --trim_tail1 ${params.tl_crop} \
              -i ${fastqs[0]} -I ${fastqs[1]} \
              -o ${sample}_paired_R1.fastq.gz -O ${sample}_paired_R2.fastq.gz \
              --unpaired1 ${sample}_unpaired.fastq.gz

        mv fastp.html fastp_${sample}.html
        """

    }
} else {
    RAW_FASTQ.trimming.set{TRIMMED_FASTQ}
}

TRIMMED_FASTQ.multiMap{it ->
    qc: it
    assembly: [it[1][0], it[1][1], it[2]]
    coverage: it[0..1]
}.set{FILTERED_FASTQ}

process PostQC {
    tag {"postQC_${sample}"}
    publishDir params.outdir+"/2-postQC", mode: "copy", pattern: "*.html"

    input:
    tuple val(sample), file(paired_fq), file(unpaired_fq) from FILTERED_FASTQ.qc

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
    mkdir paired unpaired
    mv *_paired_R*_fastqc.* paired/
    mv *_unpaired_fastqc.* unpaired/
    multiqc -n multiqc_report_paired.html paired/
    multiqc -n multiqc_report_unpaired.html unpaired/
    """
}

if (params.assembler == 'metaspades') {
    process MetaspadesAssembly {
	tag {"metaSPAdes"}
	publishDir params.outdir+"/3-assemblies", mode: "copy"

	input:
	each assembly_mode from spades_modes
	each input_mode from input_modes
	file f from FILTERED_FASTQ.assembly.collect()

	output:
	tuple(val("${assembly_mode[0]}_${input_mode}"), file("assembly_*.fasta")) into CONTIGS
    
	script:
	"""
	cat *_paired_R1.fastq.gz > paired_R1.fastq.gz
	cat *_paired_R2.fastq.gz > paired_R2.fastq.gz
	cat *_unpaired.fastq.gz > unpaired.fastq.gz

	[ "${input_mode}" = "paired-only" ] && input_mode="" || input_mode="-s unpaired.fastq.gz"
	spades.py ${assembly_mode[1]} -k 21,33,55,77 \
	    -1 paired_R1.fastq.gz -2 paired_R2.fastq.gz \$input_mode \
	    -t ${task.cpus} -m ${task.memory.getGiga()} -o spades_output

	mv spades_output/scaffolds.fasta "assembly_${assembly_mode[0]}_${input_mode}.fasta"
	"""    
    }
} else {
    process MegahitAssembly {
	tag {"megahit"}
	publishDir params.outdir+"/3-assemblies", mode: "copy"

	input:
	file f from FILTERED_FASTQ.assembly.collect()

	output:
	tuple(val('megahit'), file("assembly_*.fasta")) into CONTIGS
    
	script:
	"""
        fwd=\$(ls -m *_paired_R1.fastq.gz | tr -d ' ' | tr -d '\\n')
        rev=\$(ls -m *_paired_R2.fastq.gz | tr -d ' ' | tr -d '\\n')

	megahit -o megahit_out -1 \$fwd -2 \$rev

	mv megahit_out/final.contigs.fa assembly_megahit.fasta
	"""    
    }
    
}
    
CONTIGS.multiMap{it ->
    coverage: it
    coconet: it
    metabat2: it
    quast: it[1]
    virsorter: it
    annotation: it
}.set{ASSEMBLY}

process CoverageIndex {
    tag {"coverageIndex-${mode}"}

    input:
    tuple(val(mode), file(ctg)) from ASSEMBLY.coverage

    output:
    tuple(val(mode), file('db*')) into BWA_DB

    script:
    """
    bwa index ${ctg} -p db
    """
}

process Coverage {
    tag {"coverage-${sample}-${mode}"}
    publishDir params.outdir+"/4-coverage", mode: "copy"

    input:
    tuple(val(sample), file(paired_trimmed_fq), val(mode), file(index)) from FILTERED_FASTQ.coverage.combine(BWA_DB)

    output:
    tuple(val(mode), file('coverage*.bam')) into BAM
    
    script:
    """
    bwa mem -a -M -t 50 db ${paired_trimmed_fq} \
        | samtools sort -@ ${task.cpus} -o coverage_${mode}_${sample}.bam
    """    
}

BAM.groupTuple().multiMap{it ->
    quast: it[1]
    coconet: it
    metabat2: it}.set{COVERAGE}

process Quast {
    tag {"quast"}
    publishDir params.outdir+"/5-assemblyQC", mode: "copy"
    errorStrategy 'ignore'

    input:
    file f from ASSEMBLY.quast.collect()
    // file g from COVERAGE.quast.collect()

    output:
    file("*")

    script:
    def sp_modes = spades_modes.collect{it[0]}.join(' ')
    def in_modes = input_modes.join(' ')
    """
    # for sp_mode in $sp_modes; do
    #   for input_mode in $in_modes; do
    #     mode=\${sp_mode[0]}_\${input_mode} 
    #     samtools merge coverage_\${mode}_merged.bam coverage_\${mode}_*.bam
    #   done
    # done

    # bams=\$(ls -m *merged.bam | tr -d ' ' | tr -d '\\n')
       
    quast.py -t ${task.cpus} -o output *.fasta
    find output -type f -exec mv {} . \\;
    """    
}

process Download_db {
    tag {"virsorter-db"}

    output:
    file('virsorter-data') into virsorter_db

    script:
    """
    wget -qO- https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz | tar xz
    """    
}

process VirSorter {
    tag {"virsorter-${mode}"}
    publishDir params.outdir+"/5-assemblyQC", mode: "copy"
    container 'cyverse/virsorter'

    input:
    tuple(val(mode), file(fasta)) from ASSEMBLY.virsorter
    file(db) from virsorter_db

    output:
    file('virsorter_output/*')

    script:
    """
    wrapper_phage_contigs_sorter_iPlant.pl -f ${fasta} --db 1 \
        --wdir virsorter_output \
        --ncpu ${task.cpus} \
        --data-dir ${db}
    """
}

process CoCoNet {
    tag {"binning-coconet-${mode}"}
    publishDir params.outdir+"/6-binning", mode: "copy"
    errorStrategy 'ignore'

    input:
    tuple(val(mode), file(fasta), file(bams)) from ASSEMBLY.coconet.join(COVERAGE.coconet)

    output:
    file('coconet_bins*.csv') into COCONET_BINS

    script:
    """
    coconet --fasta ${fasta} --coverage ${bams} --output output
    mv output/bins_*.csv coconet_bins-${mode}.csv
    """
}

process Metabat2 {
    tag {"binning-metabat2-${mode}"}
    publishDir params.outdir+"/6-binning", mode: "copy"
    container "metabat/metabat"

    input:
    tuple(val(mode), file(fasta), file(bams)) from ASSEMBLY.metabat2.join(COVERAGE.metabat2)

    output:
    file('metabat2_bins*.csv') into METABAT2_BINS

    script:
    """
    runMetaBat.sh --saveCls ${fasta} ${bams}
    mv *metabat-bins*/bin metabat2_bins_${mode}.csv
    """
}

// Annotation with Prokka
process Annotation {
    tag {"annotation-${mode}"}
    publishDir params.outdir+"/7-annotation", mode: "copy"

    input:
    tuple(val(mode), file(fasta)) from ASSEMBLY.annotation

    output:
    file('prokka_out/*')

    script:
    """
    prokka ${fasta} --outdir prokka_out \
        --kingdom Viruses --metagenome \
        --mincontiglen ${params.min_annot_len} \
        --cpus ${task.cpus}
    """
}

// write summary
def summary = [:]
summary['adapters'] = params.adapters
summary['trimming tool'] = params.trimming
summary['trimming window length'] = "${params.trim_win_len} bp"
summary['trim_win window quality'] = params.trim_win_qual
summary['trimming clipping min quality'] = params.clip_qual
summary['min read length after trimming'] = "${params.min_trim_len} bp"
summary['min read average quality'] = params.min_avg_qual
summary['min contig length for annotation'] = "${params.min_annot_len} bp"
summary['Head crop'] = "${params.hd_crop} bp"
summary['Tail crop'] = "${params.tl_crop} bp"
summary['assembler'] = params.assembler

file(params.outdir).mkdir()
summary_handle = file("${params.outdir}/parameters_summary.log")
summary_handle << summary.collect { k,v -> "${k.padRight(50)}: $v" }.join("\n")
