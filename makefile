clean:
	rm -rf outputs_metagenome-assembly-pipeline work .nextflow*
test:
	nextflow main.nf --reads "test_data/*_R{1,2}.fastq.gz"
