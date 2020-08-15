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
    bedtools \
    python3-distutils \
    python3-apt

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python3 get-pip.py
RUN pip3 --version

# install numpy and pysam
WORKDIR /usr/local/bin
RUN pip3 install numpy
RUN pip3 install cython
RUN pip3 install pysam

# get script
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/anno_merger/Merge_delly_lumpy_destruct_annotation.py .
RUN cp dockerized_scripts/anno_merger/Merge_delly_lumpy_destruct_annotation.py /home/Merge_delly_lumpy_destruct_annotation.py
RUN chmod +x Merge_delly_lumpy_destruct_annotation.py
