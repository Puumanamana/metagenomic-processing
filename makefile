ifndef RELEASE
override RELEASE = latest
endif

ifndef MODULE
override MODULE = trimming
endif

clean:
	rm -rf outputs_metagenome-assembly-pipeline work .nextflow* output_test-*
test_module:
	nextflow modules/$(MODULE)/main.nf -entry test --outdir output_test-$(MODULE) -profile test
test_pipeline:
	nextflow main.nf -entry test --outdir output_test-wgs --with-singularity
test:
	make test_module MODULE=trimming
	make test_module MODULE=assembly
	make test_module MODULE=coverage
	make test_module MODULE=binning
	make test_module MODULE=annotation

container:
	cd environment \
	&& sudo docker build -f conda_dockerfile -t nakor/metagenome-assembly:$(RELEASE) . \
	&& sudo docker push nakor/metagenome-assembly:$(RELEASE) \
	&& rm -rf  ~/.singularity_images.cache/nakor-metagenome-assembly*.img
