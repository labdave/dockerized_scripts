# Base Image
FROM r-base:3.6.0

# update the OS related packages
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
    python3-pip \
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
    samtools \
    bedtools


# install pandas
WORKDIR /usr/local/bin
RUN pip3 install pandas

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("optparse", "tidyverse"), repos="http://cran.us.r-project.org")'

# get SV aggregator scripts
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/sv_merger/* .
