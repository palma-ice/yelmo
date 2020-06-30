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
    && ./configure --prefix=${LIS_PATH} --enable-omp --enable-f90 \
    && make \
    && make check \
    && make install \
    && make installcheck  \
    && rm -rf /tmp/lis-${LIS_VERSION} 

# Add user
ENV USER=glacier
RUN adduser --disabled-password --gecos '' ${USER} \
    && adduser ${USER} sudo \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER ${USER}
ENV HOME=/home/${USER}

WORKDIR /home/glacier/yelmo
ENTRYPOINT [ "bash" ]