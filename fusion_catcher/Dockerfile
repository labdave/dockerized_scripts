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
    zip \
    cmake \
    python \
    python-dev \
    python-biopython \
    gawk \
    gcc \
    g++ \
    make \
    automake \
    python-pip \
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
    bedtools \
    default-jdk \
    python-openpyxl \
    python-xlrd

# install fusion catcher
RUN wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -y -t --download

RUN python --version
RUN ls -ltr /opt/fusioncatcher/bin
RUN /opt/fusioncatcher/bin/fusioncatcher.py -h

# get fusion script
ADD https://api.github.com/repos/DevangThakkar/translocation_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/translocation_scripts.git
RUN cp translocation_scripts/fusion_catcher/fusion_catcher.sh .
RUN cp translocation_scripts/fusion_catcher/fusion_catcher.sh /home/fusion_catcher.sh
RUN chmod +x fusion_catcher.sh
