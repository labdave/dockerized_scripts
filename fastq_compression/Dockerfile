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

WORKDIR /usr/local/bin

# get spring
RUN git clone https://github.com/shubhamchandak94/Spring.git
WORKDIR /usr/local/bin/Spring
RUN mkdir build
WORKDIR /usr/local/bin/Spring/build
RUN cmake ..
RUN make
WORKDIR /usr/local/bin/
RUN cp Spring/build/spring .
RUN chmod +777 spring