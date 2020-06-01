clean:
	rm -rf outputs_metagenome-assembly-pipeline work .nextflow*
test25k:
	nextflow main.nf --reads "test_data/test_25k/*_R{1,2}.fastq.gz"
test1M:
	nextflow main.nf --reads "test_data/test_1M/*_R{1,2}.fastq.gz"
