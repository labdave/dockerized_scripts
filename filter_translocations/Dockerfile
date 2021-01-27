# Base Image
FROM r-base:3.6.0

# Metadata
LABEL base.image="filter_translocations"
LABEL description="A bioinformatics tool to filter translocation calls"
LABEL tags="Translocation Structural Variant"

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update -y && apt-get install -y \
    build-essential \
    curl \
    vim \
    less \
    wget \
    unzip \
    cmake \
    gcc-8-base \
    libmpx2 \
    libgcc-8-dev \
    libc6-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    gawk \
    python3-pip \
    bzip2 \
    git \
    autoconf \
    bsdmainutils \
    bedtools

# install Python libraries
WORKDIR /usr/local/bin
RUN pip3 install argparse
RUN pip3 install pysam

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("optparse", "tidyverse"), repos="http://cran.us.r-project.org")'

# get script
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp -r dockerized_scripts/filter_translocations/* .
RUN chmod +x filter_translocations.R

