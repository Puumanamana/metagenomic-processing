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
#                    Virsorter                    #
#-------------------------------------------------#

RUN apt-get install -y hmmer mcl
RUN wget -qO- http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz \
    | tar -xz \
    && mv mga_linux_ia64 /usr/local/bin/mga
RUN wget -qO- http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz \
    | tar -xz \
    && mv muscle3.8.31_i86linux32 /usr/local/bin/muscle
RUN wget -qO- http://github.com/bbuchfink/diamond/releases/download/v0.9.36/diamond-linux64.tar.gz \
    | tar xz \
    && mv diamond /usr/local/bin
RUN wget -qO- ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz \
    | tar -xz \
    && mv ncbi-blast-2.10.1+ /opt
RUN wget -qO- https://github.com/simroux/VirSorter/archive/v1.0.6.tar.gz \
    | tar xz \
    && cd VirSorter-1.0.6/Scripts && make clean && make && cd ../.. \
    && mv VirSorter-1.0.6 /opt \
    && ln -s /opt/VirSorter-1.0.6/wrapper_phage_contigs_sorter_iPlant.pl /usr/local/bin
RUN cpanm -l --force Bio::Perl@1.007002 Parallel::ForkManager@1.17 List::MoreUtils@0.428
RUN cpan File::Which

#-------------------------------------------------#
#    Environment variables and work directory     #
#-------------------------------------------------#

ENV PATH="/opt/FastQC:/opt/SPAdes-3.14.1-Linux/bin:${PATH}:/opt/prokka/bin:/opt/ncbi-blast-2.10.1+/bin:/opt/VirSorter-1.0.6/Scripts"

WORKDIR /workspace
COPY . /workspace

CMD /bin/bash
