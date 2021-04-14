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
    default-jdk \
    parallel

WORKDIR /usr/local/bin

# install picard
RUN ls -l
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.2/picard.jar
RUN ls -l

# get filter chrom reads script
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/filter_long_insert_exonic_reads/filter_long_insert_exonic_reads.sh .
RUN cp dockerized_scripts/filter_long_insert_exonic_reads/filter_long_insert_exonic_reads.sh /home/filter_long_insert_exonic_reads.sh
RUN chmod +x filter_long_insert_exonic_reads.sh
