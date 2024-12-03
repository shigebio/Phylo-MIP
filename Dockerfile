FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# 必要なパッケージをインストール
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

# Python 3.7 のインストール
# ubuntuのバージョンの都合上、apt-getだと取れないのでwgetしている
RUN wget https://www.python.org/ftp/python/3.7.15/Python-3.7.15.tgz && \
    tar -xzf Python-3.7.15.tgz && \
    cd Python-3.7.15 && \
    ./configure --enable-optimizations --with-openssl=/usr && \
    make && make install && \
    rm -f /usr/bin/python3 && \
    ln -s /usr/local/bin/python3.7 /usr/bin/python3 && \
    cd .. && rm -rf Python-3.7.15.tgz Python-3.7.15

# Locales の設定
RUN locale-gen en_GB.UTF-8
ENV LANG=en_GB.UTF-8
ENV LANGUAGE=en_GB:en
ENV LC_ALL=en_GB.UTF-8
ENV QT_QPA_PLATFORM=offscreen

# 環境変数設定
ENV PYTHONPATH="/usr/local/lib/python3.7/site-packages:$PYTHONPATH"

# pip のアップグレード
RUN /usr/local/bin/python3 -m ensurepip && \
    /usr/local/bin/python3 -m pip install --no-cache-dir --upgrade pip

# Python3必要なパッケージをインストール
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# PTPをクローンしてビルド
RUN git clone https://github.com/zhangjiajie/PTP /app/PTP && \
    cd /app/PTP && \
    pip3 install -r requirements.txt && \
    python3 setup.py install

# mPTPをクローンしてビルド・インストール
RUN git clone https://github.com/Pas-Kapli/mptp.git /app/mptp && \
    cd /app/mptp && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install

RUN mkdir -p /input
RUN mkdir -p /output

# ワーキングディレクトリの設定
COPY ./app /app

WORKDIR /app
