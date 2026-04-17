FROM ubuntu:22.04

LABEL maintainer="tanghaibao@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        autoconf \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        python3 \
        python3-pip \
        python3-dev \
        libjsoncpp-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir --upgrade pip \
    && python3 -m pip install --no-cache-dir boto3 awscli pandas

# Install HTSLIB / samtools
ADD install.sh /
RUN bash /install.sh

# Install splithunter, run `update_package.sh` first
ADD Splithunter /Splithunter
WORKDIR /Splithunter
RUN cd src && make -j 4 && cd ..
RUN python3 -m pip install --no-cache-dir .
