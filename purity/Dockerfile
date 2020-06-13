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
    bsdmainutils \
    bedtools

# install numpy and pysam
WORKDIR /usr/local/bin
RUN pip install numpy
RUN pip install cython
RUN pip install pysam

# get lumpy script
ADD https://api.github.com/repos/DevangThakkar/translocation_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/translocation_scripts.git
RUN cp translocation_scripts/all_merger/Merge_delly_lumpy_destruct.py .
RUN cp translocation_scripts/all_merger/Merge_delly_lumpy_destruct.py /home/Merge_delly_lumpy_destruct.py
RUN chmod +x Merge_delly_lumpy_destruct.py