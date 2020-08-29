nextflow.enable.dsl = 2

process reshape_summary {
    publishDir params.outdir
    container 'nakor/metagenome-assembly'
    
    input:
    file(summary)

    output:
    file('summary_table.csv')

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    def format_nb(n):
        if pd.isnull(n) or n<1000:
            return n
        return '{:,}k'.format(int(n/1000))

    table = pd.read_csv("${summary}", index_col=0, header=None, names=['step', 'sample', 'count'])
    
    table = (table.drop('bwa').dropna(axis=1, how='all')
        .pivot('sample', 'step', 'count')
        .reindex(columns=table.step.unique())
        .applymap(format_nb))

    table.loc['bwa'].rename(columns=['min_length', 'min_qual', 'flags'])

    table.to_csv('summary_table.csv')
    """
}

