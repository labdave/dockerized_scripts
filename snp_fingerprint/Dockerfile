# Base Image
FROM r-base:4.2.1

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update -y && apt-get install -y \
    less \
    wget \
    unzip \
    bzip2 \
    git
# install Python libraries
WORKDIR /usr/local/bin


# install R required dependencies
RUN R --vanilla -e 'install.packages(c("stringr"), repos="http://cran.us.r-project.org")'

# clone repo
ADD https://api.github.com/repos/devangthakkar/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/devangthakkar/dockerized_scripts.git

# add repo to SYSPATH
ENV PATH dockerized_scripts/snp_fingerprint:$PATH

# change the permission of the repo
RUN chmod 777 -R dockerized_scripts/snp_fingerprint
WORKDIR /usr/local/bin/dockerized_scripts/snp_fingerprint
