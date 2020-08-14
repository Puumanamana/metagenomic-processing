nextflow.enable.dsl = 2

def save_parameters() {
    def summary = [:]
    summary['Virsorter database'] = params.virsorter_db
    summary['Pfam database'] = params.pfam_db
    summary['Kraken database'] = params.kraken_db
    summary['Min length for annotation'] = params.min_annot_len
    
    file(params.outdir).mkdir()
    summary_handle = file("${params.outdir}/annotation.log")
    summary_handle << summary.collect { k,v -> "${k.padRight(50)}: $v" }.join("\n")
}

process dl_pfam_db {
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
    tuple(file(db), file(fasta))

    output:
    file('*')

    script:
    """
    viralverify.py -f ${fasta} -o ./ --hmm ${db} -t $task.cpus
    """
}

process dl_virsorter_db {
    output:
    file('virsorter-data')

    script:
    """
    wget -qO- https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz | tar xz
    """    
}

process virsorter {
    publishDir "${params.outdir}/virsorter", mode: 'copy'
    container 'cyverse/virsorter'

    input:
    tuple(file(db), file(fasta))

    output:
    file('virsorter_output/*')

    script:
    """
    wrapper_phage_contigs_sorter_iPlant.pl -f $fasta --db 1 \
        --wdir virsorter_output \
        --ncpu $task.cpus \
        --data-dir $db
    """
}

process prokka {
    publishDir "${params.outdir}/prokka", mode: 'copy'
    
    input:
    file(fasta)

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
    input:
    tuple(file(fasta), file(db))

    output:

    script:
    """
    """
}

process kraken2 {
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    
    input:
    tuple(file(fasta), file(db))

    output:

    script:
    """
    """
}

workflow annotation {
    // data is a fasta channel
    take: data
    
    main:
    if (params.virsorter_db !='null') { Channel.fromPath(params.virsorter_db) | set {vs_db} }
    else { dl_virsorter_db() | set {vs_db} }
 
    vs_db | combine(data) | virsorter

    if (params.pfam_db !='null') { Channel.fromPath(params.pfam_db) | set {pfam_db} }
    else { dl_pfam_db() | set {pfam_db} }

    pfam_db | combine(data) | viralverify

    // dl_kraken_db() | combine(data) | kraken2
    
    // data | prokka

    viralverify.out.subscribe onComplete: {save_parameters()}
}

workflow {
    Channel.fromPath(params.assembly) | annotation
}
