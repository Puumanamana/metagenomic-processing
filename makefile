ifndef RELEASE
override RELEASE = latest
endif

ifndef MODULE
override MODULE = trimming
endif

clean:
	rm -rf outputs_metagenome-assembly-pipeline work .nextflow* output_test-*
test:
	nextflow modules/$(MODULE).nf -profile lani -entry test --outdir output_test-$(MODULE)
test_all:
	make test MODULE=trimming
	make test MODULE=assembly
	make test MODULE=coverage
	make test MODULE=binning
	make test MODULE=annotation

container:
	sudo docker build -f conda_dockerfile -t nakor/metagenome-assembly:$(RELEASE) .
	sudo docker push nakor/metagenome-assembly:$(RELEASE)
	rm -rf  ~/.singularity_images.cache/nakor-metagenome-assembly*.img
