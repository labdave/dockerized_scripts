# Base Image
FROM r-base:4.0.2

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update && apt-get install -y \
    build-essential \
    #libcurl4-gnutls-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \ 
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    gawk \
    unzip \
    less \
    git \
    bedtools \
    parallel \
    samtools \
    python3 \
    python3-pip
# install Python libraries
WORKDIR /usr/local/bin


# install R required dependencies
RUN R --vanilla -e 'install.packages(c("devtools", "stringr", "plyr","doParallel", "data.table", "bedr", "optparse"), repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'devtools::install_github(repo="knausb/vcfR")'

# install pandas
RUN pip3 install pandas

# clone variantfiltering repo
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git

# add variantfiltering repo to SYSPATH
ENV PATH dockerized_scripts/variantfiltering:$PATH

# change the permission of the repo
RUN chmod 777 -R dockerized_scripts/variantfiltering
WORKDIR /usr/local/bin/dockerized_scripts/variantfiltering
