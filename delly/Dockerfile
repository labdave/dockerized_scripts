# work from latest LTS ubuntu release
FROM ubuntu:18.04

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install bcftools
RUN git clone git://github.com/samtools/htslib.git
RUN git clone git://github.com/samtools/bcftools.git
RUN ls -l
WORKDIR bcftools
 # The following is optional:
 #   autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
RUN make

# set environment
ENV BOOST_ROOT /usr

# install delly
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make STATIC=1 all \
    && make install

# copy commands to path
RUN cp /opt/delly/bin/delly /usr/bin
RUN cp bcftools /usr/bin