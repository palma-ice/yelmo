FROM continuumio/miniconda3

LABEL maintainer="thomas.goelles@gmail.com"

# Set the working directory to /home
WORKDIR /home

# languages
RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG en_US.utf8

# flags 
ENV DEBIAN_FRONTEND noninteractive 
ENV DEBCONF_NONINTERACTIVE_SEEN true
ENV LIS_PATH=/usr/lib/lis
ENV LIS_VERSION=2.0.20
ENV NETCDF_PATH=/usr
ENV FC=gfortran

#get libs  
RUN apt-get -y update && \ 
    apt-get install -y \
    git \
    gcc \
    gfortran \ 
    gdb \ 
    libnetcdf-dev \  
    libnetcdff-dev \
    less \ 
    make \
    wget \ 
    zip \
    unzip \ 
    sudo && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

# install lis
RUN echo "installing lis" \
    && cd /tmp/ \
    && wget https://www.ssisc.org/lis/dl/lis-${LIS_VERSION}.zip \
    && unzip lis-${LIS_VERSION}.zip \
    && cd /tmp/lis-${LIS_VERSION} \
    && echo $PWD \
    && echo "configure and make of lis" \
    && ./configure --prefix=${LIS_PATH} --enable-f90 \ 
    && make \
    && make check \
    && make install \
    && make installcheck  \
    && rm -rf /tmp/lis-${LIS_VERSION} 

# Python packages with conda
COPY .devcontainer/environment.yml /tmp/conda-tmp/
RUN if [ -f "/tmp/conda-tmp/environment.yml" ]; then /opt/conda/bin/conda env update -n base -f /tmp/conda-tmp/environment.yml; fi \
    && rm -rf /tmp/conda-tmp 

# runner
RUN git clone https://github.com/perrette/runner /tmp/runner  \
    && /opt/conda/bin/pip install /tmp/runner/

# Add user
ENV USER=glacier
RUN adduser --disabled-password --gecos '' ${USER} \
    && adduser ${USER} sudo \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER ${USER}
ENV HOME=/home/${USER}

RUN conda init

WORKDIR /home/glacier/yelmo