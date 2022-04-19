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

RUN git clone https://github.com/alexdobin/STAR.git
WORKDIR /usr/local/bin/STAR/source
RUN make STAR

WORKDIR /usr/local/bin

RUN mv STAR STAR_dir
# does not seem to work, will copy instead
ENV PATH "/usr/local/bin/STAR_dir/bin/Linux_x86_64/STAR:${PATH}"
RUN cp /usr/local/bin/STAR_dir/bin/Linux_x86_64/STAR .

RUN wget https://github.com/suhrig/arriba/releases/download/v2.2.1/arriba_v2.2.1.tar.gz
RUN tar -xzf arriba_v2.2.1.tar.gz
RUN cd arriba_v2.2.1 && make

RUN sed -i 's/"$THREADS"/16/g' arriba_v2.2.1/run_arriba.sh