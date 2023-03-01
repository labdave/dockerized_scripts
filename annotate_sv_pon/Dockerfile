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
    zip \
    cmake \
    python \
    python-dev \
    python-biopython \
    gawk \
    gcc \
    g++ \
    make \
    automake \
    python-pip \
    zlib1g \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev \
    libc6-dev \
    libtbb-dev \
    libtbb2 \
    bzip2 \
    tar \
    pigz \
    gzip \
    parallel \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    autoconf \
    bsdmainutils \
    bedtools \
    default-jdk \
    python-openpyxl \
    python-xlrd

# install pandas and pybedtools
WORKDIR /usr/local/bin
RUN pip install numpy
RUN pip install cython
RUN pip install pysam
RUN pip install pybedtools
RUN pip install pandas

# get script
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git
RUN cp -r dockerized_scripts/annotate_sv_pon/* .
RUN chmod +x annotate_sv_pon.py
