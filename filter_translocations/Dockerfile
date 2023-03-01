# Base Image
FROM r-base:4.0.2

# Metadata
LABEL base.image="filter_translocations"
LABEL description="A bioinformatics tool to filter translocation calls"
LABEL tags="Translocation Structural Variant"

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update && apt-get install -y \
    gawk \
    unzip \
    less \
    git \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    bedtools \
    parallel \
    samtools \
    python3 \
    python3-pip \
    libcairo2-dev \
    libxt-dev

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("optparse", "dplyr", "readr", "magrittr"), repos="http://cran.us.r-project.org")'

# get script
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git
RUN cp -r dockerized_scripts/filter_translocations/* .
RUN chmod +x filter_translocations.R

