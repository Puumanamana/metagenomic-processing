FROM python:3.7
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
	apt-transport-https \
	software-properties-common \
	cmake \
	wget \
	git \
	unzip \
	procps \
	default-jre \
	libncurses-dev \
	libbz2-dev \
	liblzma-dev \
	pkg-config \
	libfreetype6-dev \
	libpng-dev \
	libdatetime-perl \
	libxml-simple-perl \
	libdigest-md5-perl \
	bioperl

#-------------------------------------------------#
#                  Python packages                #
#-------------------------------------------------#

RUN pip3 install --upgrade pip
RUN pip3 install ipython coconet-binning multiqc

#-------------------------------------------------#
#                    Samtools                     #
#-------------------------------------------------#

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
	&& tar -xvjf samtools-1.10.tar.bz2 && cd samtools-1.10 \
	&& ./configure --prefix=/usr/local \
	&& make && make install \
	&& cd ..

#-------------------------------------------------#
#                    Bedtools                     #
#-------------------------------------------------#

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary \
	&& chmod +x bedtools.static.binary \
	&& mv bedtools.static.binary /usr/local/bin/bedtools

#-------------------------------------------------#
#                       BWA                       #
#-------------------------------------------------#

RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
	&& tar -xvjf bwa-0.7.17.tar.bz2 \
	&& cd bwa-0.7.17 && make && mv bwa /usr/local/bin \
	&& cd ..

#-------------------------------------------------#
#                  Minimap2                       #
#-------------------------------------------------#

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 \
	&& tar -xvjf minimap2-2.17.tar.bz2 \
	&& cd minimap2-2.17 && make \
	&& mv minimap2 /usr/local/bin \
	&& cd ..

#-------------------------------------------------#
#                      FastQC                     #
#-------------------------------------------------#

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
	&& unzip fastqc_v0.11.8.zip && rm fastqc_v0.11.8.zip \
	&& chmod a+x FastQC/fastqc \
	&& mv FastQC /opt

#-------------------------------------------------#
#                  Trimmomatic                    #
#-------------------------------------------------#

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
	&& unzip Trimmomatic-0.39.zip \
	&& mv Trimmomatic-0.39 /opt

#-------------------------------------------------#
#                     Fastp                       #
#-------------------------------------------------#

RUN wget https://github.com/OpenGene/fastp/archive/v0.20.1.zip \
	&& unzip v0.20.1.zip && cd fastp-0.20.1 \
	&& make && make install \
	&& cd ..

#-------------------------------------------------#
#                    Spades                       #
#-------------------------------------------------#

RUN wget https://github.com/ablab/spades/releases/download/v3.14.1/SPAdes-3.14.1-Linux.tar.gz \
	&& tar -xvzf SPAdes-3.14.1-Linux.tar.gz \
	&& mv SPAdes-3.14.1-Linux /opt

#-------------------------------------------------#
#                     Quast                       #
#-------------------------------------------------#

RUN git clone --depth 1 https://github.com/ablab/quast.git \
	&& cd quast && python3 setup.py install \
	&& cd ..

#-------------------------------------------------#
#                    Megahit                      #
#-------------------------------------------------#

RUN wget -qO- https://github.com/voutcn/megahit/archive/v1.2.9.tar.gz \
        | tar -xzf - && cd megahit-1.2.9 \
        && mkdir build && cd build \
        && cmake .. -DCMAKE_BUILD_TYPE=Release \
        && make -j4 \
        && make install

#-------------------------------------------------#
#                    Prokka                       #
#-------------------------------------------------#

RUN cpan App::cpanminus
RUN cpanm -l perl5lib Time::Piece XML::Simple Digest::MD5 Module::Build
RUN cpanm -l perl5lib Bio::AlignIO Bio::Root::Version Bio::SearchIO Bio::Seq Bio::SeqFeature::Generic Bio::SeqIO Bio::Tools::CodonTable Bio::Tools::GFF Bio::Tools::GuessSeqFormat --force
RUN git clone https://github.com/tseemann/prokka.git /opt/prokka
RUN /opt/prokka/bin/prokka --setupdb

#-------------------------------------------------#
#    Environment variables and work directory      #
#-------------------------------------------------#

ENV PATH="/opt/FastQC:/opt/SPAdes-3.14.1-Linux/bin:${PATH}:/opt/prokka/bin"

WORKDIR /workspace
COPY . /workspace

CMD /bin/bash
