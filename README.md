
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hotspots

<!-- badges: start -->

[![Travis-CI Build
Status](https://travis-ci.org/fboehm/hotspots.svg?branch=master)](https://travis-ci.org/fboehm/hotspots)

<!-- badges: end -->

The goal of hotspots is to analyze the Keller et al. (2018) expression
trait hotspots.

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

## Keller, et al. (2018)’s hotspots table

Here is Keller’s Table S1 (from the Supplemental
Documents).

| Chr | position (Mb) | number of transcripts | candidate mediator | number of genes with LOD difference \> 1.5 |
| --- | ------------- | --------------------- | ------------------ | ------------------------------------------ |
| 2   | 165.5         | 147                   | Hnf4a              | 88                                         |
| 5   | 146           | 182                   | Pdx1               | 77                                         |
| 7   | 46            | 123                   | Fam83e             | 96                                         |
| 11  | 71            | 126                   | Sat2               | 115                                        |
| 13  | 112.5         | 104                   | Il6st              | 82                                         |

From [Keller et al. 2018
GENETICS](https://www.genetics.org/content/209/1/335), [Supplementary
table 1](https://figshare.com/articles/Supplemental_Material_for_Attie_et_al_2018_in_review_/5977459)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(hotspots)
## basic example code
```

## Build the manuscript `paper.Rnw`

To render the manuscript `paper.Rnw`, into pdf, use a shell to enter the
following text from the ‘hotspots’ directory:

``` bash
docker build .
```

If you’re not using Docker, and have cloned the repository from Github,
you can use the Makefile to render `paper.Rnw`. Just type `make` in the
shell (from within the `hotspots` directory):

``` bash
make
```

## Code of Conduct

Please note that the ‘hotspots’ project is released with a [Contributor
Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.
