
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qtl2hotspots

<!-- badges: start -->

[![Travis-CI Build
Status](https://travis-ci.org/fboehm/qtl2hotspots.svg?branch=master)](https://travis-ci.org/fboehm/qtl2hotspots)

<!-- badges: end -->

The goal of `qtl2hotspots` is to analyze the Keller et al. (2018)
expression trait hotspots with both mediation methods and pleiotropy
tests.

## Installation

Install from Dockerhub with this command:

``` bash
docker pull fjboehm/qtl2hotspots
```

Alternatively, if you wish to install the package directly from Github,
use `devtools`.

``` r
install.packages("devtools")
devtools::install_github("fboehm/qtl2hotspots")
```

## Keller, et al. (2018)’s hotspots table

Here is Keller’s Table S1 (from the Supplemental Documents).

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
library(tidyverse)
#> ── Attaching packages ──────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
#> ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
#> ✔ tibble  2.1.3     ✔ dplyr   0.8.3
#> ✔ tidyr   0.8.3     ✔ stringr 1.4.0
#> ✔ readr   1.3.1     ✔ forcats 0.4.0
#> ── Conflicts ─────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
```

``` r
map <- readRDS(here::here("analysis", "data", "derived_data", "map.rds"))
hot_chr <- c(2, 5, 7, 11, 13)
start_index <- c(3742, 3498, 927, 1868, 2814) # from s1 values in condor submit files
scan_length <- c(247, 269, 242, 217, 226) # in number of markers, from condor submit files
stop_index <- start_index + scan_length - 1 # index of the last marker in the interval
start_position <- c(map$`2`[start_index[1]], 
                    map$`5`[start_index[2]], 
                    map$`7`[start_index[3]],
                    map$`11`[start_index[4]],
                    map$`13`[start_index[5]]
)
stop_position <- c(map$`2`[stop_index[1]], 
                    map$`5`[stop_index[2]], 
                    map$`7`[stop_index[3]],
                    map$`11`[stop_index[4]],
                    map$`13`[stop_index[5]]
)
tibble::tibble(hot_chr, start_index, scan_length, stop_index, start_position, stop_position) %>%
  knitr::kable()
```

| hot\_chr | start\_index | scan\_length | stop\_index | start\_position | stop\_position |
| -------: | -----------: | -----------: | ----------: | --------------: | -------------: |
|        2 |         3742 |          247 |        3988 |       162.79779 |       168.3768 |
|        5 |         3498 |          269 |        3766 |       142.49079 |       148.8150 |
|        7 |          927 |          242 |        1168 |        40.75175 |        50.3110 |
|       11 |         1868 |          217 |        2084 |        67.60676 |        76.7000 |
|       13 |         2814 |          226 |        3039 |       109.57503 |       115.7887 |

## Build the manuscript `paper.Rnw`

To render the manuscript `paper.Rnw`, into pdf, use a shell to enter the
following text from the `qtl2hotspots` directory:

``` bash
docker build .
```

If you’re not using Docker, and have cloned the repository from Github,
you can use the Makefile to render `paper.Rnw`. Just type `make` in the
shell (from within the `qtl2hotspots` directory):

``` bash
make
```

## Code of Conduct

Please note that the `qtl2hotspots` project is released with a
[Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By
contributing to this project, you agree to abide by its terms.

## Status

Sept 10, 2019: Might need to redo 2d scans with bigger scan regions.
