# work from amcpherson's docker image
FROM amcpherson/destruct:0.4.17

# reset entrypoint
ENTRYPOINT ["/usr/bin/env"]

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
    zip

# install samtools
ENV samtools_version 1.9
WORKDIR /usr/local/bin/
RUN curl -SL https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    > /usr/local/bin/samtools-${samtools_version}.tar.bz2
RUN tar -xjf /usr/local/bin/samtools-${samtools_version}.tar.bz2 -C /usr/local/bin/
RUN cd /usr/local/bin/samtools-${samtools_version}/ && ./configure
RUN cd /usr/local/bin/samtools-${samtools_version}/ && make
RUN cd /usr/local/bin/samtools-${samtools_version}/ && make install

# change reference
ADD https://api.github.com/repos/labdave/dockerized_scripts/git/refs/heads/ version.json
RUN git clone https://github.com/labdave/dockerized_scripts.git

# location 1
RUN mkdir -p temp/
RUN unzip /opt/conda/pkgs/destruct-0.4.17-py37_1/lib/python3.7/site-packages/destruct-0.4.17-py3.7.egg -d temp/
RUN rm /opt/conda/pkgs/destruct-0.4.17-py37_1/lib/python3.7/site-packages/destruct-0.4.17-py3.7.egg
RUN cp dockerized_scripts/destruct/hg38_chr_map.tsv temp/destruct/data/
RUN cp dockerized_scripts/destruct/download.py temp/destruct/utils/
RUN cp dockerized_scripts/destruct/defaultconfig.py temp/destruct/

RUN rm temp/destruct/utils/__pycache__/download.cpython-37.pyc
RUN python -m compileall temp/destruct/utils/download.py
RUN rm temp/destruct/__pycache__/defaultconfig.cpython-37.pyc
RUN python -m compileall temp/destruct/defaultconfig.py

WORKDIR temp/
RUN zip -r destruct-0.4.17-py3.7.egg destruct EGG-INFO
WORKDIR /usr/local/bin/
RUN cp temp/destruct-0.4.17-py3.7.egg /opt/conda/pkgs/destruct-0.4.17-py37_1/lib/python3.7/site-packages/
WORKDIR /usr/local/bin/

# location 2
RUN mkdir -p temp1/
RUN unzip /opt/conda/lib/python3.7/site-packages/destruct-0.4.17-py3.7.egg -d temp1/
RUN rm /opt/conda/lib/python3.7/site-packages/destruct-0.4.17-py3.7.egg
RUN cp dockerized_scripts/destruct/hg38_chr_map.tsv temp1/destruct/data/
RUN cp dockerized_scripts/destruct/download.py temp1/destruct/utils/
RUN cp dockerized_scripts/destruct/defaultconfig.py temp1/destruct/

RUN rm temp1/destruct/utils/__pycache__/download.cpython-37.pyc
RUN python -m compileall temp1/destruct/utils/download.py
RUN rm temp1/destruct/__pycache__/defaultconfig.cpython-37.pyc
RUN python -m compileall temp1/destruct/defaultconfig.py

WORKDIR temp1/
RUN zip -r destruct-0.4.17-py3.7.egg destruct EGG-INFO
WORKDIR /usr/local/bin/
RUN cp temp/destruct-0.4.17-py3.7.egg /opt/conda/lib/python3.7/site-packages/
WORKDIR /usr/local/bin/

# get ref data
RUN destruct create_ref_data destruct_ref/