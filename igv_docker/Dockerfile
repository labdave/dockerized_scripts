# work from latest LTS ubuntu release
FROM openjdk:8-jdk-slim

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
    xvfb \
    libxext6 \
    libxrender-dev \
    libfontconfig1 \
    libxrender1 \
    libxtst6 \
    x11-utils

# get igv
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.3/IGV_2.3.81.zip && unzip IGV_2.3.81.zip
ADD https://api.github.com/repos/DevangThakkar/translocation_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/igv_docker/test.script .

# this is crucial - would simply not work without this - I don't pretend to know why
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]