# Base Image
FROM r-base:4.0.2

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
    parallel \
WORKDIR /usr/local/bin

# clone genotyping repo
ADD https://api.github.com/repos/devangthakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/devangthakkar/dockerized_scripts.git

# add genotyping repo to SYSPATH
ENV PATH dockerized_scripts/genotyping:$PATH

# change the permission of the repo
RUN chmod 777 -R dockerized_scripts/genotyping
WORKDIR /usr/local/bin/dockerized_scripts/genotyping
