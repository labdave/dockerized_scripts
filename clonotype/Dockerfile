################################################################################
##################### Set Inital Image to work from ############################
FROM ghcr.io/milaboratory/mixcr/mixcr:4.0.0-68-develop

################################################################################
##################### Add Container Labels #####################################
LABEL "Description"="Software toolkit for T/B cell repertoire"

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
##################### download clonotype calling scripts #######################

ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git

# download BCR script
RUN mv dockerized_scripts/clonotype/BCR.sh .

# download TCR script
RUN mv dockerized_scripts/clonotype/TCR.sh .

# change entrypoint
ENTRYPOINT ["/bin/bash", "-c"]
