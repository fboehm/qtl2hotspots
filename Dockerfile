# get the base image, the rocker/verse has R, RStudio and pandoc
FROM fjboehm/verse-with-pandoc-crossref

# required - see https://github.com/thomasWeise/docker-pandoc/blob/98ad0ab53308701b93fa8ba82c80546495fc7828/image/Dockerfile
MAINTAINER Frederick J. Boehm <frederick.boehm@gmail.com>

COPY . /home/rstudio/qtl2hotspots

  # build this compendium package

RUN cd /home/rstudio/qtl2hotspots; make -f Makefile-docker.makefile
