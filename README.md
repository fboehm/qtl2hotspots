
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hotspots

<!-- badges: start -->

[![Travis-CI Build
Status](https://travis-ci.org/fboehm/hotspots.svg?branch=master)](https://travis-ci.org/fboehm/hotspots)

<!-- badges: end -->

The goal of hotspots is to …

## Installation

Install from Dockerhub with this command:

``` bash
docker pull fjboehm/hotspots
```

Alternatively, if you wish to install the package directly from Github,
use `devtools`.

``` r
install.packages("devtools")
devtools::install_github("fboehm/hotspots")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(hotspots)
## basic example code
```

## Build the manuscript `paper.Rmd`

To render the manuscript `paper.Rmd`, into html, tex, and pdf, use a
shell to enter the following text from the ‘hotspots’ directory:

``` bash
docker build .
```

If you’re not using Docker, and have cloned the repository from Github,
you can use the Makefile to render `paper.Rmd`. Just type `make` in the
shell:

``` bash
make
```

## Code of Conduct

Please note that the ‘hotspots’ project is released with a [Contributor
Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.
