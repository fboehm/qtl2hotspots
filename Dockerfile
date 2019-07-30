# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/verse:3.5.3

# required - see https://github.com/thomasWeise/docker-pandoc/blob/98ad0ab53308701b93fa8ba82c80546495fc7828/image/Dockerfile
MAINTAINER Frederick J. Boehm <frederick.boehm@gmail.com>

COPY . /home/rstudio/qtl2hotspots

# Install R packages
RUN install2.r --error \
    broman \
    here \
    heatmaply \
    dendextend \
    plotly \
    svglite \
    bookdown \
    webshot \
    knitr \
    purrr \
    dplyr \
    htmlwidgets \
    devtools \
    git2r \
    tidyverse


RUN apt-get clean &&\
    apt-get update &&\
    apt-get autoclean -y &&\
    apt-get autoremove -y &&\
    apt-get update &&\
# install tzdata, which we will need down the line
    apt-get install -f -y tzdata &&\
# We install cabal, a Haskell package manager, because we want the newest
# pandoc and filters which we can only get from there.
# We also install zlib1g, as we will need it later on.
# We install librsvg2 in order to make svg -> pdf conversation possible.
    apt-get install -f -y haskell-platform \
                          cabal-install \
                          libudunits2-dev && \
# get the newest list of packages
    cabal update &&\
# install the dependencies of the packages we want
    cabal install --dependencies-only \
                  pandoc \
                  pandoc-citeproc \
                  pandoc-crossref &&\
# install the packages we want
    cabal install pandoc \
                  pandoc-citeproc \
                  pandoc-crossref &&\
# clean up all temporary files
    apt-get clean &&\
    apt-get autoclean -y &&\
    apt-get autoremove -y &&\
    apt-get clean

# we remember the path to pandoc in a special variable
ENV PANDOC_DIR=/root/.cabal/bin/

# add pandoc to the path
ENV PATH=${PATH}:${PANDOC_DIR}

# go into the repo directory
RUN . /etc/environment \
  # Install linux dependencies here
  # build this compendium package
  && cd /home/rstudio/qtl2hotspots; make
