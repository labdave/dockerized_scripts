# work from r release
FROM r-base:4.0.2

MAINTAINER DaveLab <lab.dave@gmail.com>

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
    python3 \
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
    rna-star \
    samtools \
    bwa \
    bsdmainutils

# install numpy and pysam
WORKDIR /usr/local/bin
RUN pip3 install numpy
RUN pip3 install cython
RUN pip3 install pysam

#install R packages
RUN R --vanilla -e 'install.packages(c("tibble"), repos="http://cran.us.r-project.org")'

# add proper davrelab docker GitHub resources
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git
#add script for aggregating msi
RUN cp dockerized_scripts/msi/Aggregate_MSI.R /home/
