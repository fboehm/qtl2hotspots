# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/verse:3.6.1

# required - see https://github.com/thomasWeise/docker-pandoc/blob/98ad0ab53308701b93fa8ba82c80546495fc7828/image/Dockerfile
MAINTAINER Frederick J. Boehm <frederick.boehm@gmail.com>

COPY . /home/rstudio/qtl2hotspots

# Install packages
RUN DEBIAN_FRONTEND=noninteractive apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yqq \
    libssl-dev \
    libssh2-1-dev \
    libffi-dev \
    zlib1g-dev \
    python-cffi \
    python-dev \
    python-pip \
    build-essential \
    cmake \
    gcc \
    pkg-config \
    git \
    libhttp-parser-dev \
    python-setuptools \
    wget

RUN wget https://github.com/libgit2/libgit2/archive/v0.28.0.tar.gz && \
tar xzf v0.28.0.tar.gz && \
cd libgit2-0.28.0/ && \
cmake . && \
make && \
make install

RUN ldconfig



# Install R packages
RUN install2.r --error \
    BiocManager

RUN Rscript -e 'BiocManager::install("S4Vectors")'

RUN install2.r --error \
    broman \
    here \
    dendextend \
    plotly \
    git2r \
    iheatmapr \
    pheatmap


# build this compendium package
RUN cd /home/rstudio/qtl2hotspots/analysis/paper/Rmd; make
