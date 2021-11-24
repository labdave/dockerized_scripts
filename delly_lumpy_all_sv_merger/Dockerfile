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
ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git
RUN cp dockerized_scripts/delly_lumpy_all_sv_merger/merge_delly_lumpy_all_SVs.py .
RUN cp dockerized_scripts/delly_lumpy_all_sv_merger/merge_delly_lumpy_all_SVs.py /home/merge_delly_lumpy_all_SVs.py
