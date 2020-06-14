ifndef RELEASE
override RELEASE = latest
endif

clean:
	rm -rf outputs_metagenome-assembly-pipeline work .nextflow*
test25k:
	nextflow main.nf --reads "test_data/test_25k/*_R{1,2}.fastq.gz" -profile lani
test1M:
	nextflow main.nf --reads "test_data/test_1M/*_R{1,2}.fastq.gz" -profile lani
container:
	sudo docker build -f dockerfile -t nakor/metagenome-assembly:$(RELEASE) .
	sudo docker push nakor/metagenome-assembly:$(RELEASE)
	rm -rf  ~/.singularity_images.cache/nakor-metagenome-assembly*.img
