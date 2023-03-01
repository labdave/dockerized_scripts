# work from latest LTS ubuntu release
FROM python:3.6

# run update and install necessary tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    curl \
    vim \
    less \
    wget \
    unzip \
    zip \
    cmake \
    gawk \
    gcc \
    g++ \
    make \
    automake \
    zlib1g \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev \
    libc6-dev \
    libtbb-dev \
    libtbb2 \
    bzip2 \
    tar \
    pigz \
    gzip \
    parallel \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    autoconf \
    bsdmainutils \
    bedtools

RUN which pip3

# install numpy and pandas
WORKDIR /usr/local/bin
RUN pip3 install numpy
RUN pip3 install pandas

# get script
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git
RUN cp dockerized_scripts/anno_merger/merge_D_L_with_DW.py .
RUN cp dockerized_scripts/anno_merger/merge_D_L_with_DW.py /home/merge_D_L_with_DW.py
RUN cp dockerized_scripts/anno_merger/multisample_merge_D_L_with_DW.py .
RUN cp dockerized_scripts/anno_merger/multisample_merge_D_L_with_DW.py /home/multisample_merge_D_L_with_DW.py
