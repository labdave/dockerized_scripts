# work from latest LTS ubuntu release
FROM openjdk:18-jdk-slim

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
    python3 \
    gawk \
    python3-pip \
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
    bedtools \
    samtools \
    xvfb \
    libxext6 \
    libxrender-dev \
    libfontconfig1 \
    libxrender1 \
    libxtst6 \
    x11-utils

ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git
RUN cp -r dockerized_scripts/igv_docker/* .

# get igv
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.3/IGV_2.3.81.zip && unzip IGV_2.3.81.zip
RUN cp dockerized_scripts/igv_docker/prefs.properties IGV_2.3.81/

RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.14/IGV_2.14.1.zip && unzip IGV_2.14.1.zip
RUN cp dockerized_scripts/igv_docker/prefs.properties IGV_2.14.1/

#add pil
RUN pip3 install pillow

# this is crucial - would simply not work without this - I don't pretend to know why
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]

# Copy session file
COPY hg38_session.xml /data/hg38_session.xml
COPY Twist_8MB_panel_with_ERCCs.maskPAR.bed /data/Twist_8MB_panel_with_ERCCs.maskPAR.bed
COPY Twist_exome.maskPAR.bed /data/Twist_exome.maskPAR.bed

# work from latest LTS ubuntu release
#FROM ubuntu:18.04

# run update and install necessary tools
#RUN apt-get update -y && apt-get install -y \
#    build-essential \
#    libnss-sss \
#    curl \
#    vim \
#    less \
#    wget \
#    unzip \
#    cmake \
#    python \
#    gawk \
#    python-pip \
#    zlib1g-dev \
#    libncurses5-dev \
#    libncursesw5-dev \
#    libnss-sss \
#    libbz2-dev \
#    liblzma-dev \
#    bzip2 \
#    libcurl4-openssl-dev \
#    libssl-dev \
#    git \
#    autoconf \
#    bsdmainutils \
#    bedtools \
#    samtools \
#    xvfb \
#    libxext6 \
#    libxrender-dev \
#    libfontconfig1 \
#    libxrender1 \
#    libxtst6 \
#    x11-utils

# get igv
#RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_Linux_2.8.6.zip && unzip IGV_Linux_2.8.6.zip
#ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
#RUN git clone https://github.com/labdave/dockerized_scripts.git
#RUN cp dockerized_scripts/igv_docker/filter_rows.py .
#RUN cp dockerized_scripts/igv_docker/test.script .
#RUN cp dockerized_scripts/igv_docker/igv_script_creator_lanie.py .
#RUN cp dockerized_scripts/igv_docker/prefs.properties IGV_Linux_2.8.6/

# this is crucial - would simply not work without this - I don't pretend to know why
#ENV TINI_VERSION v0.19.0
#ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
#RUN chmod +x /tini
#ENTRYPOINT ["/tini", "--"]
