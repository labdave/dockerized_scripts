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
    python3-pip \
    libcairo2-dev \
    libxt-dev
# install Python libraries
WORKDIR /usr/local/bin


# install R required dependencies
RUN R --vanilla -e 'install.packages(c("BiocManager", "Cairo"), repos="http://cran.us.r-project.org")'
RUN R --vanilla -e 'BiocManager::install("preprocessCore")'
RUN R --vanilla -e 'BiocManager::install("ComplexHeatmap")'
RUN R --vanilla -e 'install.packages(c("FARDEEP", "ggplot2", "colorspace", "RColorBrewer", "circlize", "tidyverse"), repos="http://cran.us.r-project.org")'

# install pandas
RUN pip3 install pandas

# clone fardeep repo
ADD https://api.github.com/repos/devangthakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/devangthakkar/dockerized_scripts.git

# add fardeep repo to SYSPATH
ENV PATH dockerized_scripts/fardeep:$PATH

# change the permission of the repo
RUN chmod 777 -R dockerized_scripts/fardeep
WORKDIR /usr/local/bin/dockerized_scripts/fardeep
