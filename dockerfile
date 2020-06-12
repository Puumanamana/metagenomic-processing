FROM python:3.7
MAINTAINER Cedric Arisdakessian <carisdak@hawaii.edu>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y apt-transport-https software-properties-common wget unzip nano emacs procps default-jre

#-------------------------------------------------#
#                  Python packages                #
#-------------------------------------------------#

# RUN apt-get update && apt-get install -y python3 python3-pip
# RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1
RUN pip3 install --upgrade pip
RUN pip3 install ipython coconet-binning multiqc

#-------------------------------------------------#
#                    Samtools                     #
#-------------------------------------------------#

RUN apt-get install -y libncurses-dev libbz2-dev liblzma-dev

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

RUN apt-get install -y pkg-config libfreetype6-dev libpng-dev git
RUN git clone --depth 1 https://github.com/ablab/quast.git \
	&& cd quast && python3 setup.py install \
	&& cd ..

#-------------------------------------------------#
#                    Megahit                      #
#-------------------------------------------------#

RUN apt-get install -y cmake

RUN wget -qO- https://github.com/voutcn/megahit/archive/v1.2.9.tar.gz \
        | tar -xzf - && cd megahit-1.2.9 \
        && mkdir build && cd build \
        && cmake .. -DCMAKE_BUILD_TYPE=Release \
        && make -j4 \
        && make install

#-------------------------------------------------#
#    Environment variables and work directory      #
#-------------------------------------------------#

ENV PATH="/opt/FastQC:/opt/SPAdes-3.14.1-Linux/bin:${PATH}"

WORKDIR /workspace
COPY . /workspace

CMD /bin/bash
