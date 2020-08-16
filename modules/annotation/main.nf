nextflow.enable.dsl = 2

def save_parameters() {
    def summary = [:]
    summary['Virsorter database'] = params.virsorter2_db
    summary['Pfam database'] = params.pfam_db
    summary['Kraken database'] = params.kraken_db
    summary['Min length for annotation'] = params.min_annot_len
    
    file(params.outdir).mkdir()
    summary_handle = file("${params.outdir}/annotation.log")
    summary_handle << summary.collect { k,v -> "${k.padRight(50)}: $v" }.join("\n")
}

process dl_pfam_db {
    publishDir "$HOME/db", mode: 'copy'
    
    output:
    file('Pfam-A.hmm')

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz 
    gunzip Pfam-A.hmm.gz
    """
}

process viralverify {
    publishDir "${params.outdir}/viralverify", mode: 'copy'
    
    input:
    path fasta
    path db

    output:
    file('*')

    script:
    """
    viralverify.py -f ${fasta} -o ./ --hmm ${db} -t $task.cpus
    """
}

process dl_virsorter2_db {
    publishDir "$HOME/db", mode: 'copy'

    output:
    path 'virsorter2_db'

    script:
    """
    virsorter setup -d virsorter2_db -j $task.cpus
    """    
}

process virsorter2 {
    publishDir "${params.outdir}/virsorter", mode: 'copy'

    input:
    path fasta
    path db

    output:
    file('virsorter_output')

    script:
    """
    #!/usr/bin/env bash

    virsorter run -w virsorter_output -i $fasta -j $task.cpus -d $db
    """
}

process prokka {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/prokka", mode: 'copy'
    
    input:
    path fasta

    output:
    file('prokka_out/*')

    script:
    """
    prokka $fasta --outdir prokka_out \
        --kingdom Viruses --metagenome \
        --mincontiglen ${params.min_annot_len} \
        --centre X --compliant \
        --cpus $task.cpus
    """
}

process dl_kraken2_db {
    errorStrategy 'ignore'
    
    input:
    path fasta
    path db

    output:

    script:
    """
    """
}

process kraken2 {
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    path fasta
    path db

    output:

    script:
    """
    """
}

workflow annotation {
    take: fasta
    
    main:
    if (params.virsorter2_db !='null') { Channel.fromPath(params.virsorter2_db) | set {vs2_db} }
    else { dl_virsorter2_db() | set {vs2_db} }
 
    virsorter2(fasta, vs2_db)

    if (params.pfam_db !='null') { Channel.fromPath(params.pfam_db) | set {pfam_db} }
    else { dl_pfam_db() | set {pfam_db} }

    viralverify(fasta, pfam_db)

    // dl_kraken_db() | combine(data) | kraken2
    
    prokka(fasta)
    viralverify.out.subscribe onComplete: {save_parameters()}
}

workflow {
   annotation(params.assembly)
}

workflow test {
    annotation("$baseDir/../../test_data/assembly.fasta")
}
