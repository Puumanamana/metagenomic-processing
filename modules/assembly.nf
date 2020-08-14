nextflow.enable.dsl = 2

def save_parameters() {
    def summary = [:]
    summary['Assembler'] = params.assembler
    summary['Co-assembly'] = !params.split_assembly
    
    file(params.outdir).mkdir()
    summary_handle = file("${params.outdir}/assembler.log")
    summary_handle << summary.collect { k,v -> "${k.padRight(50)}: $v" }.join("\n")
}

process spades {
    publishDir "${params.outdir}/${params.assembler}", mode: 'copy'
    
    input:
	tuple(val(sample), file(fastq))

	output:
	tuple(val(sample), file("assembly_*.fasta"))
    
	script:
	"""
    #!/usr/bin/env bash

    inputs="-1 \$(ls *_R1.fastq.gz) -2 \$(ls *_R2.fastq.gz)"

    if compgen -G "*_unpaired.fastq.gz" > /dev/null; then
        inputs+=" -s \$(ls *_unpaired.fastq.gz)"
    fi

	${params.assembler}.py -k 21,33,55,77 \$inputs -t $task.cpus -m ${task.memory.getGiga()} -o ${params.assembler}_output

	mv ${params.assembler}_output/scaffolds.fasta "assembly_${params.assembler}_${sample}.fasta"
	"""    
    }

process megahit {
    publishDir "${params.outdir}/megahit", mode: 'copy'
    
    input:
    tuple(val(sample), file(fastq))

	output:
	tuple(val(sample), file("assembly_*.fasta"))
    
	script:
	"""
    #!/usr/bin/env bash

    inputs="-1 \$(ls *_R1.fastq.gz) -2 \$(ls *_R2.fastq.gz)"
    if compgen -G "*_unpaired.fastq.gz" > /dev/null; then
        inputs+=" -r \$(ls *_unpaired.fastq.gz)"
    fi

	megahit -o megahit_out \$inputs

	mv megahit_out/final.contigs.fa assembly_megahit_${sample}.fasta
	"""    
}

process quast {
    publishDir "${params.outdir}/quast", mode: 'copy'
    errorStrategy 'ignore'

    input:
    file(fasta)

    output:
    file("*")

    script:
    """
    quast.py -t $task.cpus -o output *.fasta
    find output -type f -exec mv {} . \\;
    """    
}


workflow assembly {
    take: data

    main:
    if(params.split_assembly) {
        reads = data
    } else {
        data | transpose | collectFile(sort: true) {
            it -> ["all_" + (it[1].getSimpleName() =~/_[^_]*/)[-1][1..-1] + ".fastq.gz", it[1]]
        } | map{['pool', it]} | groupTuple | set{reads}
    }
    
    if (params.assembler == 'megahit') {reads | megahit | set{assembly}}
    else {reads | spades | set{assembly}}
    assembly | collect | quast

    assembly.subscribe onComplete: {save_parameters()}
    
    emit: assembly
}

workflow {
    Channel.fromFilePairs(params.reads) | assembly
}

workflow test {
    Channel.fromFilePairs("$baseDir/../test_data/*_R{1,2}.fastq.gz") | assembly
}
