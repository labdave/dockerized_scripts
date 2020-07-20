# work from latest LTS ubuntu release
FROM ubuntu:18.04

# run update and install necessary tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    curl \
    vim \
    less \
    wget \
    unzip \
    cmake \
    python \
    gawk \
    python-pip \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    autoconf \
    bsdmainutils \
    bedtools \
    samtools \
    default-jre \
    default-jdk

# install numpy and pysam
WORKDIR /usr/local/bin
RUN pip install numpy
RUN pip install cython
RUN pip install pysam

# install picard
RUN ls -l
RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.2/picard.jar
RUN ls -l

# get filter chrom reads script
ADD https://api.github.com/repos/DevangThakkar/translocation_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/translocation_scripts.git
RUN cp translocation_scripts/filter_chrom_reads/filter_chrom_reads.sh .
RUN cp translocation_scripts/filter_chrom_reads/filter_chrom_reads.sh /home/filter_chrom_reads.sh
RUN chmod +x filter_chrom_reads.sh
