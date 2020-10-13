# work from latest LTS ubuntu release
FROM ubuntu:18.04

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
    python \
    gawk \
    python-pip \
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
    bsdmainutils

# install numpy and pysam
WORKDIR /usr/local/bin
RUN pip install pandas

# get cnv aggregator scripts
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/cnv_merger/get_gene_cn.py .
RUN cp dockerized_scripts/cnv_merger/get_gene_cn.py /home/get_gene_cn.py
RUN chmod +x get_gene_cn.py
RUN cp dockerized_scripts/cnv_merger/get_cyto_cn.py .
RUN cp dockerized_scripts/cnv_merger/get_cyto_cn.py /home/get_cyto_cn.py
RUN chmod +x get_cyto_cn.py
