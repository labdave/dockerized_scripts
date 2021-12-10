# Base Image
FROM r-base:4.0.2

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

# install pandas
WORKDIR /usr/local/bin
RUN pip3 install pandas

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("optparse", "dplyr", "readr", "purrr", "stringr"), repos="http://cran.us.r-project.org")'

# get SV aggregator scripts
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/sv_merger/* .
