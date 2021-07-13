################################################################################
##################### Set Inital Image to work from ############################
FROM ubuntu:16.04

################################################################################
##################### Add Container Labels #####################################
LABEL "Description"="Software toolkit for T/B cell repertoire"

################################################################################
############## set environment variables #######################################

ENV version 3.0.3

################################################################################
##################### Install System Dependencies ##############################
RUN apt-get update -y && apt-get install -y \
    git \
    build-essential \
    openjdk-8-jre \
    wget \
    unzip \
    samtools \
    bedtools

################################################################################
##################### Install MIXCR    #########################################

# set inital WORKDIR
WORKDIR /mixcr/

# download the source
RUN wget https://github.com/milaboratory/mixcr/releases/download/v${version}/mixcr-${version}.zip

# unzip the source
RUN unzip mixcr-${version}.zip

################################################################################
##################### download clonotype calling scripts #######################

ADD https://api.github.com/repos/DevangThakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/DevangThakkar/dockerized_scripts.git

# download BCR script
RUN mv dockerized_scripts/clonotype/BCR.sh .

# download TCR script
RUN mv dockerized_scripts/clonotype/TCR.sh .

################################################################################
###################### set environment path    #################################

# add mixcr executable to path
ENV PATH="/mixcr/mixcr-${version}/:${PATH}"
