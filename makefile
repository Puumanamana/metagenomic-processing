ifndef RELEASE
override RELEASE = latest
endif

ifndef MODULE
override MODULE = trimming
endif

ifndef PROFILE
override PROFILE = singularity
endif

clean:
	rm -rf outputs_metagenome-assembly-pipeline work .nextflow* output_test-*
test_module:
	nextflow modules/$(MODULE)/main.nf -entry test --outdir output_test-$(MODULE) -profile $(PROFILE)
test_pipeline:
	nextflow main.nf -entry test --outdir output_test-wgs -profile $(PROFILE)
test:
	make test_module MODULE=trimming PROFILE=$(PROFILE)
	make test_module MODULE=assembly PROFILE=$(PROFILE)
	make test_module MODULE=coverage PROFILE=$(PROFILE)
	make test_module MODULE=binning PROFILE=$(PROFILE)
	make test_module MODULE=annotation PROFILE=$(PROFILE)

container:
	cd environment \
	&& sudo docker build -f conda_dockerfile -t nakor/metagenome-assembly:$(RELEASE) . \
	&& sudo docker push nakor/metagenome-assembly:$(RELEASE) \
	&& rm -rf  ~/.singularity_images.cache/nakor-metagenome-assembly*.img
