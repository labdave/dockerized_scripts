# work from latest LTS ubuntu release
FROM openjdk:11-jdk-slim

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
RUN git clone https://github.com/stevekm/IGV-snapshot-automator.git
WORKDIR IGV-snapshot-automator/bin
RUN bash get_IGV.sh
WORKDIR /IGV-snapshot-automator