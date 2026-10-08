# SPDX-FileCopyrightText: 2024-2026 <actual copyright holder(s)>
# SPDX-License-Identifier: GPL-3.0-only

FROM ubuntu@sha256:8feb4d8ca5354def3d8fce243717141ce31e2c428701f6682bd2fafe15388214

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONPATH="/usr/local/lib/python3.7/site-packages"

COPY LICENSE /app/LICENSE
COPY THIRD_PARTY_NOTICES.md /app/THIRD_PARTY_NOTICES.md

# Install required packages
RUN apt-get update && apt-get install -y \
    lsb-release=11.1.0ubuntu2 \
    build-essential=12.8ubuntu1.1 \
    git=1:2.25.1-1ubuntu3.14 \
    cmake=3.16.3-1ubuntu1.20.04.1 \
    unzip=6.0-25ubuntu1.2 \
    wget=1.20.3-1ubuntu2.1 \
    vsearch=2.14.1-3build1 \
    mafft=7.453-1 \
    fasttree=2.1.11-1 \
    libcurl4-openssl-dev=7.68.0-1ubuntu2.25 \
    libxml2-dev=2.9.10+dfsg-5ubuntu0.20.04.10 \
    software-properties-common=0.99.9.12 \
    libreadline-dev=8.0-4 \
    libpcre++-dev=0.9.5-6.1build1 \
    libblas-dev=3.9.0-1build1 \
    liblapack-dev=3.9.0-1build1 \
    libatlas-base-dev=3.10.3-8ubuntu7 \
    gfortran=4:9.3.0-1ubuntu2 \
    locales=2.31-0ubuntu9.18 \
    autoconf=2.69-11.1 \
    automake=1:1.16.1-4ubuntu6 \
    flex=2.6.4-6.2 \
    bison=2:3.5.1+dfsg-1 \
    g++=4:9.3.0-1ubuntu2 \
    xvfb=2:1.20.13-1ubuntu1~20.04.20 \
    libtool=2.4.6-14 \
    qt5-default=5.12.8+dfsg-0ubuntu2.1 \
    libqt5x11extras5=5.12.8-0ubuntu1 \
    libxkbcommon-x11-0=0.10.0-1 \
    libxcb-xinerama0=1.14-2 \
    libgl1-mesa-dev=21.2.6-0ubuntu0.1~20.04.2 \
    libgl1-mesa-glx=21.2.6-0ubuntu0.1~20.04.2 \
    liblzma-dev=5.2.4-1ubuntu1.1 \
    libbz2-dev=1.0.8-2 \
    python3-pyqt5=5.14.1+dfsg-3build1 \
    libhdf5-dev=1.10.4+repack-11ubuntu1 \
    zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.5 \
    libffi-dev=3.3-4 \
    libsqlite3-dev=3.31.1-4ubuntu0.7 \
    libxslt1-dev=1.1.34-4ubuntu0.20.04.3 \
    python3-dev=3.8.2-0ubuntu2 \
    libssl-dev=1.1.1f-1ubuntu2.24 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

# Installing Python 3.7
ARG PYTHON_SOURCE_SHA256=cf2993798ae8430f3af3a00d96d9fdf320719f4042f039380dca79967c25e436
RUN wget https://www.python.org/ftp/python/3.7.15/Python-3.7.15.tgz && \
    echo "${PYTHON_SOURCE_SHA256}  Python-3.7.15.tgz" | sha256sum -c - && \
    tar -xzf Python-3.7.15.tgz && \
    cd Python-3.7.15 && \
    ./configure --enable-optimizations --with-openssl=/usr && \
    make && make install && \
    rm -f /usr/bin/python3 && \
    ln -s /usr/local/bin/python3.7 /usr/bin/python3 && \
    cd .. && rm -rf Python-3.7.15.tgz Python-3.7.15

# Setting Locales
RUN locale-gen en_GB.UTF-8
ENV LANG=en_GB.UTF-8
ENV LANGUAGE=en_GB:en
ENV LC_ALL=en_GB.UTF-8
ENV QT_QPA_PLATFORM=offscreen

# Setting environment variables
ENV PYTHONPATH="/usr/local/lib/python3.7/site-packages:$PYTHONPATH"

# Upgrading pip
RUN /usr/local/bin/python3 -m ensurepip && \
    /usr/local/bin/python3 -m pip install --no-cache-dir pip==24.0 setuptools==47.1.0 wheel==0.42.0

# Install the packages in requirements.txt
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Clone and build bPTP
ARG PTP_COMMIT=f8d1da0888bf2ce84d4f8ce160915f51439c94ec
RUN git clone https://github.com/zhangjiajie/PTP /app/PTP && \
    cd /app/PTP && \
    git checkout ${PTP_COMMIT} && \
    pip3 install -r requirements.txt && \
    python3 setup.py install

# Clone and build mPTP
ARG MPTP_COMMIT=1f98d29aac4c6ccd5a2737412891b450c71a8480
RUN git clone https://github.com/Pas-Kapli/mptp.git /app/mptp && \
    cd /app/mptp && \
    git checkout ${MPTP_COMMIT} && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install

# Copy application files
COPY ./app /app

# Copy entrypoint script
COPY ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

USER root

# Make scripts executable in container
RUN echo '#!/bin/bash\npython3 /app/Phylo-MIP.py "$@"' > /usr/local/bin/phylo-mip && \
    echo '#!/bin/bash\npython3 /app/merge_data.py "$@"' > /usr/local/bin/merge_data && \
    chmod +x /usr/local/bin/phylo-mip && \
    chmod +x /usr/local/bin/merge_data

ENV PATH="/usr/local/bin:${PATH}"

WORKDIR /app

# Ensure we always execute Python scripts with python3
ENTRYPOINT ["/entrypoint.sh"]
