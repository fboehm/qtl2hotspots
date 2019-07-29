# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/verse:3.5.3

# required
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




# go into the repo directory
RUN . /etc/environment \
  \
  # Install linux depedendencies here
  && sudo apt-get update \
  && sudo apt-get install libudunits2-dev -y \
  && sudo apt-get install cabal-install \
  && cabal update && cabal install pandoc-crossref \
  # build this compendium package
  && cd /home/rstudio/qtl2hotspots; make
