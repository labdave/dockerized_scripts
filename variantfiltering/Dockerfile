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
    bedtools \
    parallel \
    samtools \
    python3 \
    python3-pip
# install Python libraries
WORKDIR /usr/local/bin


# install R required dependencies
RUN R --vanilla -e 'install.packages(c("devtools", "stringr", "plyr","doParallel", "data.table", "bedr"), repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'devtools::install_github(repo="knausb/vcfR")'

# install pandas
RUN pip3 install pandas

# clone variantfiltering repo
ADD https://api.github.com/repos/devangthakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/devangthakkar/dockerized_scripts.git

# add variantfiltering repo to SYSPATH
ENV PATH dockerized_scripts/variantfiltering:$PATH

# change the permission of the repo
RUN chmod 777 -R dockerized_scripts/variantfiltering
WORKDIR /usr/local/bin/dockerized_scripts/variantfiltering
