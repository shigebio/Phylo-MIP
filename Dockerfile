FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONPATH="/usr/local/lib/python3.7/site-packages"

# Install required packages
RUN apt-get update && apt-get install -y \
    lsb-release \
    build-essential \
    git \
    cmake \
    unzip \
    wget \
    vsearch \
    mafft \
    fasttree \
    libcurl4-openssl-dev \
    libxml2-dev \
    software-properties-common \
    libreadline-dev \
    libpcre++-dev \
    libblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    gfortran \
    locales \
    autoconf \
    automake \
    flex \
    bison \
    g++ \
    xvfb \
    libtool \
    qt5-default \
    libqt5x11extras5 \
    libxkbcommon-x11-0 \
    libxcb-xinerama0 \
    libgl1-mesa-dev \
    libgl1-mesa-glx \
    liblzma-dev \
    libbz2-dev \
    python3-pyqt5 \
    libhdf5-dev \
    zlib1g-dev \
    libffi-dev \
    libsqlite3-dev \
    libxslt1-dev \
    python3-dev \
    libssl-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

# Installing Python 3.7
RUN wget https://www.python.org/ftp/python/3.7.15/Python-3.7.15.tgz && \
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
    /usr/local/bin/python3 -m pip install --no-cache-dir --upgrade pip

# Install the packages in requirements.txt
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Clone and build bPTP
RUN git clone https://github.com/zhangjiajie/PTP /app/PTP && \
    cd /app/PTP && \
    pip3 install -r requirements.txt && \
    python3 setup.py install

# Clone and build mPTP
RUN git clone https://github.com/Pas-Kapli/mptp.git /app/mptp && \
    cd /app/mptp && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install

# Copy application files
COPY ./app /app

# Copy entrypoint script
COPY ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Make scripts executable in container
RUN echo '#!/bin/bash\npython3 /app/MICUM.py "$@"' > /usr/local/bin/micum && \
    echo '#!/bin/bash\npython3 /app/merge_data.py "$@"' > /usr/local/bin/merge_data && \
    chmod +x /usr/local/bin/micum && \
    chmod +x /usr/local/bin/merge_data

ENV PATH="/usr/local/bin:${PATH}"

WORKDIR /app

# Ensure we always execute Python scripts with python3
ENTRYPOINT ["/entrypoint.sh"]