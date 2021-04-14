# Base Image
FROM r-base:3.6.0

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
    gawk \
    bzip2 \
    git \
    autoconf \
    bsdmainutils

# install pandas
WORKDIR /usr/local/bin
RUN pip install pandas

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("optparse", "tidyverse"), repos="http://cran.us.r-project.org")'

# get SV aggregator scripts
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/sv_merger/* .
