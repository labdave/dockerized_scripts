# work from latest LTS ubuntu release
FROM ubuntu:18.04

# set the environment variables
ENV lumpy_version 0.3.0
ENV samblaster_version 0.1.24
ENV sambamba_version 0.6.9
ENV samtools_version 1.9

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
    bsdmainutils

# install numpy and pysam
WORKDIR /usr/local/bin
RUN pip install numpy
RUN pip install cython
RUN pip install pysam

# get lumpy script
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git
RUN cp dockerized_scripts/lumpy_merger/Merge_sample_level_Lumpy.py .
RUN cp dockerized_scripts/lumpy_merger/Merge_sample_level_Lumpy.py /home/Merge_sample_level_Lumpy.py
RUN chmod +x Merge_sample_level_Lumpy.py
