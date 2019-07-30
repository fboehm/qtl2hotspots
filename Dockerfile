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
    webshot \
    htmlwidgets \
    git2r \
    bookdown


RUN apt-get clean &&\
    apt-get update &&\
    apt-get autoclean -y &&\
    apt-get autoremove -y &&\
    apt-get update &&\
# install tzdata, which we will need down the line
    apt-get install -f -y tzdata &&\
# We install cabal, a Haskell package manager, because we want the newest
# pandoc and filters which we can only get from there.
    apt-get install -f -y libudunits2-dev
#### install haskell stuff -- from https://github.com/freebroccolo/docker-haskell/blob/7fd359b8dab3bf543832eb1ff34e1a46eef262a7/8.6/Dockerfile
RUN apt-get update && \
    apt-get install -y --no-install-recommends gnupg ca-certificates dirmngr curl git && \
    echo 'deb http://downloads.haskell.org/debian stretch main' > /etc/apt/sources.list.d/ghc.list && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys BA3CBA3FFE22B574 && \
    apt-get update && \
    apt-get install -y --no-install-recommends ghc-8.6.5 cabal-install-2.4 \
        zlib1g-dev libtinfo-dev libsqlite3-dev g++ netbase xz-utils make && \
    curl -fSL https://github.com/commercialhaskell/stack/releases/download/v1.9.3/stack-1.9.3-linux-x86_64.tar.gz -o stack.tar.gz && \
    curl -fSL https://github.com/commercialhaskell/stack/releases/download/v1.9.3/stack-1.9.3-linux-x86_64.tar.gz.asc -o stack.tar.gz.asc && \
    export GNUPGHOME="$(mktemp -d)" && \
    gpg --batch --keyserver ha.pool.sks-keyservers.net --recv-keys C5705533DA4F78D8664B5DC0575159689BEFB442 && \
    gpg --batch --verify stack.tar.gz.asc stack.tar.gz && \
    tar -xf stack.tar.gz -C /usr/local/bin --strip-components=1 && \
    /usr/local/bin/stack config set system-ghc --global true && \
    /usr/local/bin/stack config set install-ghc --global false && \
    rm -rf "$GNUPGHOME" /var/lib/apt/lists/* /stack.tar.gz.asc /stack.tar.gz

ENV PATH /root/.cabal/bin:/root/.local/bin:/opt/cabal/2.4/bin:/opt/ghc/8.6.5/bin:$PATH




# get the newest list of packages
RUN cabal update
# install the dependencies of the packages we want
RUN cabal install --dependencies-only \
                  pandoc \
                  pandoc-citeproc \
                  pandoc-crossref
# install the packages we want
RUN cabal install pandoc \
                  pandoc-citeproc \
                  pandoc-crossref \
# we remember the path to pandoc in a special variable
ENV PANDOC_DIR=/root/.cabal/bin/
# add pandoc to the path
ENV PATH=${PATH}:${PANDOC_DIR}
# go into the repo directory
RUN . /etc/environment \
  # Install linux dependencies here
  # build this compendium package
  && cd /home/rstudio/qtl2hotspots; make
