Expression QTL hotspot dissection with mediation analysis and pleiotropy
testing
================
Frederick J. Boehm, Brian S. Yandell, and Karl W. Broman
7/27/2019

# Introduction

Genetics studies in model organisms like mice can identify genomic
regions that affect complex quantitative traits, such as systolic blood
pressure, body weight, and, more recently, biomolecular traits such as
gene expression levels and protein concentrations (Sax
[1923](#ref-sax1923association); Soller, Brody, and Genizi
[1976](#ref-soller1976power); Lander and Botstein
[1989](#ref-lander1989mapping); Broman and Sen
[2009](#ref-broman2009guide); Jansen
[2007](#ref-jansen2007quantitative); Chick et al.
[2016](#ref-chick2016defining)). These genomic regions are called
“quantitative trait loci” or “QTL”. A genome-wide QTL “scan” reveals
associations between genotypes and phenotypes by considering each
position along the genome, one at a time, as a candidate QTL for the
trait of interest. A region with strong evidence of association with a
complex quantitative trait, then, defines a QTL (for that trait).
Because nearby markers have correlated genotypes, a QTL in a two-parent
cross often spans multiple megabases in length and may contain more than
a hundred genes. Identification of a causal gene (for a given complex
trait) from among those genes near the QTL is challenging and may
require costly and time-consuming experiments. The growing need for
greater QTL mapping resolution fueled development over the last two
decades of model organism multiparental populations for high-resolution
QTL mapping (Koning and McIntyre [2017](#ref-de2017back); Churchill et
al. [2004](#ref-churchill2004collaborative); Svenson et al.
[2012](#ref-svenson2012high); Huang et al.
[2012](#ref-huang2012multiparent), [2011](#ref-huang2011analysis);
Shivakumar et al. [2018](#ref-shivakumar2018soybean); Kover et al.
[2009](#ref-kover2009multiparent); Tisne et al.
[2017](#ref-tisne2017identification); Stanley et al.
[2017](#ref-stanley2017genetic)).

Recent systems genetics studies have sought to identify QTL for
thousands of biomolecular traits (Keller et al.
[2018](#ref-keller2018genetic); Chick et al.
[2016](#ref-chick2016defining)).

## INTRO IDEAS

  - why is it useful to study expression hotspots?
  - who else has studied hotspots? And what did they do?

<!-- end list -->

``` r
library(broman)
library(here)
```

# Methods

Keller et al. ([2018](#ref-keller2018genetic)) identified five
expression trait hotspots in their study of pancreatic islet gene
expression in 378 Diversity Outbred mice. They defined a hotspot to be a
4-Mb region that affects at least 100 nonlocal expression traits.

For each hotspot, we identified a set of local expression traits and a
set of nonlocal expression traits. All traits demonstrated LOD peaks
that met the genome-wide significance threshold, 7.18 (Keller et al.
[2018](#ref-keller2018genetic)). The nonlocal expression traits arose
from genes on other chromosomes, while the local traits reflected
expression levels for genes on the same chromosome as the hotspot and
within 4 Mb of the hotspot center. We ignored transcripts for genes on
the hotspot chromosome but not within 4 Mb of the hotspot center.

## Pleiotropy testing

### Two-dimensional, two-QTL scans

<span id="eq:mvlmm">\[
vec(Y) = X vec(B) + vec(G) + vec(E)
\qquad(1)\]</span>

For each hotspot, we tested all pairs of traits, regardless of whether
the traits were local or nonlocal. We performed a two-dimensional,
two-QTL scan over a genomic region that included the entirety of the
hotspot. We fitted bivariate linear mixed effects models (Equations
[1](#eq:mvlmm)-[3](#eq:e)) at each ordered pair of markers in a
two-dimensional grid (F. J. Boehm et al. [2019](#ref-boehm2019testing)).

<span id="eq:g">\[
G \sim MN_{n x 2}(0, K, V_g)
\qquad(2)\]</span>

<span id="eq:e">\[
E \sim MN_{n x 2}(0, I_n, V_e)
\qquad(3)\]</span>

and \(G\) and \(E\) are independent.

\[
\hat B = (X^T\hat \Sigma^{-1}X)^{-1}X^T\hat\Sigma^{-1}Y
\]

\[
\hat \Sigma = \hat V_g \otimes \hat K + \hat V_e \otimes I_n
\]

\[
l_{10}(\hat B, \hat \Sigma) = - \frac{1}{2}\log_{10}\left((Y - X\hat B)^T\hat \Sigma^{-1}(Y - X\hat B)\right) - n\log_{10}(2\pi) - \log_{10}|\hat \Sigma|
\]

For \(n\) mice, \(\Sigma\) is a \(2n\) by \(2n\) covariance matrix,
while \(V_g\) and \(V_e\) are \(2\) by \(2\) covariance matrices.
Estimates \(\hat V_g\) and \(\hat V_e\) are obtained via an
expectation-maximization algorithm (Zhou and Stephens
[2014](#ref-zhou2014efficient); Boehm [2018](#ref-gemma2)).

### Calculating likelihood ratio test statistics

From the two-dimensional scans, we calculated pleiotropy likelihood
ratio test statistics for each pair of traits (Jiang and Zeng
[1995](#ref-jiang1995multiple); F. J. Boehm et al.
[2019](#ref-boehm2019testing); F. J. Boehm
[2019](#ref-boehm2019thesis)). Because of the required computing time
for bootstrap analyses, we didn’t obtain p-values for the pleiotropy
tests. Instead, we worked directly with the test statistic values.

<span id="eq:lod">\[
LOD = \log_{10} \left(\frac{\max_{\lambda_1, \lambda_2}L(B, \Sigma, \lambda_1, \lambda_2)}{\max_{\lambda}L(B, \Sigma, \lambda, \lambda)} \right)
\qquad(4)\]</span>

\[
ll(\lambda_1, \lambda_2)
\]

## Mediation analysis

A mediation analysis in systems genetics begins with an identified
QTL-expression trait association. In our case, we choose a QTL that
affects a nonlocal trait.

For each hotspot and each local trait - nonlocal trait pair, we also
performed a mediation analysis to assess the extent which the local
trait mediates the QTL effect on the nonlocal trait. Each mediation
analysis involved fitting four linear models (Equations\~, , , ) (Chick
et al. [2016](#ref-chick2016defining)).

<span id="eq:model1">\[
Y = \alpha 1 + WC + E
\qquad(5)\]</span>

<span id="eq:model2">\[
Y = XB + WC + E
\qquad(6)\]</span>

<span id="eq:model3">\[
Y = \alpha 1 + WC + M\beta + E
\qquad(7)\]</span>

<span id="eq:model4">\[
Y = XB + WC + M\beta + E
\qquad(8)\]</span>

In all of our mediation analyses, the nonlocal transcripts are the
targets, while the local transcripts serve as candidate mediators of the
QTL-target relationship. We only considered as potential mediators those
local traits with strong univariate LODs. In Equations\~ through , \(X\)
is a matrix of founder allele dosages, \(B\) is a matrix of founder
allele effects, \(W\) is a matrix of covariates, \(C\) is a matrix of
covariate effects, \(M\) is a vector of values for a single candidate
mediator, \(\beta\) is the effect of the mediator on the complex trait,
1 is a vector with all entries set to 1, \(\alpha\) is the mean trait
value, and \(E\) is a vector of normally distributed random errors.

## Visualizations

We created scatter plots with the `ggplot2` package (Wickham
[2016](#ref-ggplot2)) for the R statistical computing environment (R
Core Team [2018](#ref-r)). We initially plotted all local-nonlocal pairs
for each hotspot. We highlighted in blue those points that correspond to
pairs involving the putative mediator per (Keller et al.
[2018](#ref-keller2018genetic)).

### Heatmaps

We created heatmaps to examine patterns among pleiotropy test results in
each hotspot. Our analysis used both local and nonlocal traits. We
represented the pleiotropy test statistics as a symmetric matrix, where
each row was one expression trait and each column was an expression
trait. Excepting the cells on the diagonal, each cell, then, contained
the pleiotropy test statistic from the analysis involving the row’s
expression trait and the column’s expression trait. To create a heatmap,
we first performed two hierarchical clusterings of rows. One
hierarchical clustering involved only the rows and columns containing
local traits, while the other contained only rows and columns containing
nonlocal traits. We calculated the Euclidean distances between all row
pairs (in the respective submatrices) before clustering. Hierarchical
clustering was performed using the complete linkage method. Initially,
each row is assigned to its own cluster. An iterative algorithm then
combines the two most similar clusters at each step. Distances between
clusters were computed with the Lance-Williams dissimilarity update
formula for complete linkage (Lance and Williams
[1967](#ref-lance1967general)). Lastly, we arranged the rows according
to optimal leaf ordering (Hahsler, Hornik, and Buchta
[2008](#ref-hahsler2008getting); Galili [2015](#ref-dendextend)). We
then combined the local and nonlocal traits into a single heatmap for
each hotspot. To achieve coloring contrasts in the heatmaps, we set the
color spectrum limits at 0 (dark blue) and 5 (yellow). Pleiotropy test
statistics above 5, too, are represented by yellow. We annotated rows by
a collection of columns on the right-hand side of the heatmap. This
binary matrix consists of columns that summarize mediation results. Each
column represents a local trait. If a row corresponds to a nonlocal
trait that is mediated by the local trait (*i.e.*, with mediation LOD
difference greater than 1.5), then the cell is colored dark green;
otherwise, the cell is light green. We didn’t examine the possibility
that a local trait might mediate the relationship between a QTL and
another local trait; thus, all local trait rows are colored light green.
The last annotation column and the annotation row above the heatmap both
designate local traits with dark blue and nonlocal traits with light
blue.

# Results

We generated scatter plots of pleiotropy test statistics vs. mediation
LOD difference values for each of the five hotspots (Figures). Each
point represents a pairing of a local trait with a nonlocal trait. The
mediation LOD difference value is that observed when the local trait is
considered as a candidate mediator of the QTL-nonlocal trait
association.

In all five scatter plots, we observe few points in the upper right
quadrant of the figure. That is, few pairs have high pleiotropy test
statistics and high mediation LOD difference values.

``` r
library(qtl2hotspots)
pleio_threshold <- 4
ult <- 15 # univariate lod threshold
```

``` r
library(tidyverse)
options(tibble.print_max = 50, tibble.print_min = 20)
```

``` r
library(heatmaply)
```

``` r
if (!file.exists(here("analysis", "data", "derived_data", "expr.rds"))){
  download.file("https://datadryad.org/bitstream/handle/10255/dryad.166774/Attie_DO378_eQTL_viewer_v1.Rdata?sequence=2",
              destfile = here("analysis", "data", "raw_data", "Attie_DO378_eQTL_viewer_v1.Rdata")
              )
  load(here("analysis", "data", "raw_data", "Attie_DO378_eQTL_viewer_v1.Rdata"))
  saveRDS(dataset.islet.rnaseq$lod.peaks, here("analysis", "data", "derived_data", "lod_peaks.rds"))
  saveRDS(dataset.islet.rnaseq$annots, here("analysis", "data", "derived_data", "annots.rds"))
  saveRDS(dataset.islet.rnaseq$expr, here("analysis", "data", "derived_data", "expr.rds"))
  saveRDS(K, here("analysis", "data", "derived_data", "kinship.rds"))
  saveRDS(map, here("analysis", "data", "derived_data", "map.rds"))
  saveRDS(genoprobs, here("analysis", "data", "derived_data", "genoprobs.rds"))
  saveRDS(dataset.islet.rnaseq$covar, here("analysis", "data", "derived_data", "covar.rds"))
}
```

``` r
lod_peaks <- readRDS(here("analysis", "data", "derived_data", "lod_peaks.rds"))
annots <- readRDS(here("analysis", "data", "derived_data", "annots.rds"))
expr <- readRDS(here("analysis", "data", "derived_data", "expr.rds"))
```

``` r
hot_chr <- 2
hot_mid <- 165.5
keller_mediator <- "Hnf4a"
```

# Chr 2

``` r
hot_lower <- hot_mid - 2
hot_upper <- hot_mid + 2
fp <- paste0("Chr", hot_chr, "-")
knitr::opts_chunk$set(fig.path = here("analysis", "figures", fp))
```

``` r
pleio_mixed_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio.rds")))
pleio_local_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio_local.rds")))
pleio_nonlocal_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio_nonlocal.rds")))
med_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_mediations.rds")))
```

``` r
local_annot <- lod_peaks %>%
  filter(chrom == hot_chr, pos <= hot_upper, pos >= hot_lower) %>%
  left_join(y = annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= ult) %>%
  mutate(local_tx = chr == hot_chr) %>%
  filter(local_tx)
```

``` r
nonlocal_annot <- lod_peaks %>%
  filter(chrom == hot_chr, pos <= hot_upper, pos >= hot_lower) %>%
  inner_join(y = annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= 7.18) %>%
  mutate(local_tx = chr == hot_chr) %>%
  filter(!local_tx) %>%
  filter(!is.na(hotspot))
```

``` r
local_annot2 <- local_annot %>%
  filter(middle <= hot_upper, middle >= hot_lower) %>% 
  arrange(desc(lod))
```

``` r
hot_local <- expr[, colnames(expr) %in% local_annot2$annot.id]
hot_nonlocal <- expr[, colnames(expr) %in% nonlocal_annot$annot.id]
nnonlocal <- nrow(nonlocal_annot)
nlocal <- nrow(local_annot2)
```

``` r
med_tib <- med_list %>%
  purrr::map(bind_rows) %>%
  bind_rows() %>%
  as_tibble() %>%
  rename(lod_without_med = lod_no_med, lod_with_med = lod_med) %>%
  mutate(lod_diff = lod_without_med - lod_with_med, lod_diff_proportion = lod_diff / lod_without_med) %>%
  left_join(local_annot2, by = c("local_gene_id" = "annot.id")) %>%
  left_join(nonlocal_annot, by = c("nonlocal_gene_id" = "annot.id")) %>%
  dplyr::select(symbol.x, symbol.y, lod_without_med, lod_with_med, lod_diff, lod_diff_proportion, local_gene_id, nonlocal_gene_id) %>%
  rename(local_symbol = symbol.x, nonlocal_symbol = symbol.y) %>%
  filter(!is.na(local_symbol))
```

``` r
hot_annot <- local_annot2 %>%
  bind_rows(nonlocal_annot) %>%
  dplyr::select(annot.id, symbol)
pleio_tib <- pleio_mixed_list %>%
  purrr::map(.f = function(x) rename(x, gene1_id = local_gene_id, gene2_id = nonlocal_gene_id)) %>%
  c(pleio_local_list, pleio_nonlocal_list) %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(hot_annot, by = c("gene1_id" = "annot.id")) %>%
  rename(gene1_symbol = symbol) %>%
  inner_join(hot_annot, by = c("gene2_id" = "annot.id")) %>%
  rename(gene2_symbol = symbol) %>%
  dplyr::select(gene1_symbol, gene2_symbol, pleiotropy_lod)
```

``` r
pm <- pleio_tib %>%
  tibble_to_matrix(symmetric = TRUE)
```

``` r
p <- pleio_mixed_list %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(med_tib, by = c("local_gene_id", "nonlocal_gene_id")) %>%
  dplyr::select(local_symbol, nonlocal_symbol, pleiotropy_lod, lod_diff, lod_diff_proportion, lod_without_med, lod_with_med, filename, nonlocal_gene_id, local_gene_id) %>%
  ggplot() + geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = local_symbol == keller_mediator, nonlocal = nonlocal_symbol, local = local_symbol), size = 0.05) + theme(legend.position = "none") + xlab("Pleiotropy LOD") + ylab("Mediation LOD difference")
```

    ## Warning: Ignoring unknown aesthetics: nonlocal, local

``` r
fnsvg <- paste0("Chr", hot_chr, "_scatter.svg")
fnpdf <- paste0("Chr", hot_chr, "_scatter.pdf")
ggsave(plot = p, filename = here("analysis", "figures", fnsvg))
```

    ## Saving 7 x 5 in image

``` r
ggsave(plot = p, filename = here("analysis", "figures", fnpdf))
```

    ## Saving 7 x 5 in image

``` r
plotly::ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr2-unnamed-chunk-19-1.png)<!-- -->

``` r
foo <- pleio_mixed_list %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(med_tib, by = c("local_gene_id", "nonlocal_gene_id")) %>%
  dplyr::select(local_symbol, nonlocal_symbol, pleiotropy_lod, lod_diff, lod_diff_proportion, lod_without_med, lod_with_med, filename, nonlocal_gene_id, local_gene_id)
p <- foo %>%
  ggplot() +
  geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = nonlocal_symbol), size = 0.1) +
  facet_wrap(. ~ local_symbol, nrow = ceiling(nlocal / 2), ncol = 2) +
  broman::karl_theme(strip.placement = "inside") +
  theme(legend.position = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
) +
  xlab("Pleiotropy LOD") +
  ylab("Mediation LOD difference") +
  geom_text(aes(x = max(foo$pleiotropy_lod), y = 1, label = local_symbol), inherit.aes=FALSE)

fnsvg <- paste0("Chr", hot_chr, "_scatter-panel.svg")
fnpdf <- paste0("Chr", hot_chr, "_scatter-panel.pdf")
ggsave(p, filename = here("analysis", "figures", fnsvg), height = ceiling(nlocal / 2))
```

    ## Saving 7 x 16 in image

``` r
ggsave(p, filename = here("analysis", "figures", fnpdf), height = ceiling(nlocal / 2))
```

    ## Saving 7 x 16 in image

``` r
plotly::ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr2-unnamed-chunk-20-1.png)<!-- -->

``` r
# put in multiple figures, if there are many local traits
n_figs <- 1 + nlocal %/% 18 # 18 is the number of panels per page in pdf
n_last <- nlocal %% 18 # number of panels on last figure
max_pleio_lod <- max(foo$pleiotropy_lod)
max_med_lod_difference <- max(foo$lod_diff)
min_med_lod_difference <- min(foo$lod_diff)

for (i in 1:n_figs){
  n_panels <- 18 * (i < n_figs) + n_last * (i == n_figs)
  foo %>%
    filter(local_symbol %in% local_annot2$symbol[((i - 1) * 18 + 1): (18 * i)]) %>%
    ggplot() +
      geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = nonlocal_symbol), size = 0.1) +
      facet_wrap(. ~ local_symbol, nrow = ceiling(n_panels / 2), ncol = 2) +
      broman::karl_theme() +
      theme(legend.position = "none") +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
      xlab("Pleiotropy LOD") + xlim(c(0, max_pleio_lod)) + ylim(c(min_med_lod_difference, max_med_lod_difference)) +
      ylab("Mediation LOD difference") +
      geom_text(size = 2, aes(x = max_pleio_lod / 2, y = 2 * max_med_lod_difference / 3, label = local_symbol), inherit.aes=FALSE)
  fnsvg <- paste0("Chr", hot_chr, "_scatter-panel", i, ".svg")
  fnpdf <- paste0("Chr", hot_chr, "_scatter-panel", i, ".pdf")
  ggsave(filename = here("analysis", "figures", fnsvg), height = ceiling(n_panels / 2))
  ggsave(filename = here("analysis", "figures", fnpdf), height = ceiling(n_panels / 2))
}
```

``` r
scols1 <- 1L + (rownames(pm) %in% local_annot2$symbol)
scols2 <- rownames(pm) %in% (med_tib %>% filter(local_symbol == keller_mediator, lod_diff > 1.5))$nonlocal_symbol
```

``` r
# order the nonlocal traits first
## indicator for being a nonlocal trait
indic_nonlocal <- rownames(pm) %in% nonlocal_annot$symbol
# define ord_nonlocal
ord_nonlocal <- pm[indic_nonlocal, indic_nonlocal] %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(pm[indic_nonlocal, indic_nonlocal])) %>%
  order.dendrogram()
rn_nonlocal <- rownames(pm[indic_nonlocal, indic_nonlocal][ord_nonlocal, ord_nonlocal])
# define ord_local
ord_local<- pm[!indic_nonlocal, !indic_nonlocal] %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(pm[!indic_nonlocal, !indic_nonlocal])) %>%
  order.dendrogram()
rn_local <- rownames(pm[!indic_nonlocal, !indic_nonlocal][ord_local, ord_local])
rn <- c(rn_nonlocal, rn_local)
```

``` r
scols_extra <- purrr::map(.x = local_annot2$symbol,
                          .f = function(x) {rownames(pm) %in% (med_tib %>% filter(local_symbol == x, lod_diff > 1.5))$nonlocal_symbol}
) %>%
  bind_cols() %>%
  (function(x){colnames(x) <- local_annot2$symbol; x}) %>%
  mutate(local = scols1)
```

``` r
n_med_tr <- 4 # arbitrary threshold for number of mediated traits per local
# only those local traits with n_med_tr or more mediated nonlocal traits 
# appear as annotation columns in the heatmap
n_med_per_local <- med_tib %>%
  group_by(local_symbol) %>%
  filter(lod_diff > 1.5) %>%
  count() %>%
  filter(n >= n_med_tr)
n_med_per_local %>%
  arrange(desc(n)) %>%
  knitr::kable()
```

| local\_symbol |  n |
| :------------ | -: |
| Hnf4a         | 90 |
| Slc12a5       | 39 |
| Slc35c2       | 38 |
| Pabpc1l       | 24 |
| Ctsa          | 21 |
| Slc2a10       | 21 |
| Tomm34        | 15 |
| Gm14291       | 11 |
| Eya2          | 10 |
| Snx21         | 10 |
| 1700025C18Rik |  8 |
| Neurl2        |  8 |
| Stk4          |  5 |
| Gm11460       |  4 |
| Serinc3       |  4 |

``` r
scols_pre <- scols_extra %>%
  dplyr::select(n_med_per_local$local_symbol)
# get ordering for columns of scols_pre
col_order <- scols_pre %>%
  t() %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(t(scols_pre))) %>%
  order.dendrogram()
scols <- scols_pre[ , col_order] %>%
  mutate(local = scols1)
fnhtml <- here("analysis", "figures", paste0("Chr", hot_chr, "hm.html"))
```

``` r
ord <- match(rn, rownames(pm))

(p <- pm[ord, ord] %>%
  apply(FUN = trunc2, MARGIN = 2, threshold = 5) %>%
  heatmaply(Rowv = FALSE,
            Colv = FALSE,
            cexRow = 0.2,
            cexCol = 0.2,
            symm = TRUE,
            row_side_colors = scols[ord, ],
            col_side_colors = tibble(local = scols1[ord]),
            limits = c(0, 5)
            ) %>%
  hide_guides()
)
```

    ## Warning in heatmaply.heatmapr(hm, colors = colors, limits = limits,
    ## scale_fill_gradient_fun = scale_fill_gradient_fun, : The hover text for
    ## col_side_colors is currently not implemented (due to an issue in plotly).
    ## We hope this would get resolved in future releases.

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr2-unnamed-chunk-26-1.png)<!-- -->

``` r
fnpdf <- here("analysis", "figures", paste0("Chr", hot_chr, "hm.pdf"))
htmlwidgets::saveWidget(widget = p, 
                        file = fnhtml
                        )
```

``` r
webshot::webshot(url = fnhtml,
                 file = fnpdf,
                 cliprect = "viewport"
                 )
```

``` r
bar <- foo %>%
  filter(local_symbol == keller_mediator)

p <- tibble(nonlocal_symbol = rownames(pm[indic_nonlocal, indic_nonlocal])[ord_nonlocal], index = 1:length(nonlocal_symbol)) %>%
  inner_join(bar, by = "nonlocal_symbol") %>% 
  left_join(nonlocal_annot, by = c("nonlocal_symbol" = "symbol")) %>%
  ggplot() + geom_point(aes(y = pos, x = index, color = lod_diff > 1.5, nonlocal = nonlocal_symbol, lod_diff = lod_diff)) 
```

    ## Warning: Ignoring unknown aesthetics: nonlocal, lod_diff

``` r
ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr2-unnamed-chunk-28-1.png)<!-- -->

## Chr 2 hotspot results

The scatter plot of pleiotropy test statistics and mediation LOD
differences reveals that few local-nonlocal trait pairs have both a high
value of the pleiotropy test statistic and a high value of the mediation
LOD difference. Many points have low values of both, while some points
demonstrate a high value of only mediation LOD difference or pleiotropy
test statistic.

Many points have negative values of mediation LOD difference. These
correspond to pairs that, when conditioning on the putative mediator,
the LOD for the nonlocal trait increased.

The points are color-coded to identify those pairs that involve the
local trait *Hnf4a*, which Keller et al.
([2018](#ref-keller2018genetic)) identified as a key mediator of traits
at this hotspot. Blue points involve *Hnf4a*, while red points do not.
We see that many blue point have large mediation LOD differences and
small pleiotropy test statistics. In fact, most of the points with
mediation LOD difference above 15 involve *Hnf4a* and have small
pleiotropy test statistics.

An unanticipated observation is that some points have mediation LOD
differences of nearly 10 and pleiotropy test statistics of nearly 10.
With such large pleiotropy test statistics, one has reasonable evidence
each point representing two separate QTL. It seems biologically unlikely
that a nonlocal trait and a mediating local trait would not share a
pleiotropic QTL.

In looking closer at these six points, we find that five of the six have
*Tmem189* as the local trait. This raises questions about the *Tmem189*
trait. One known issue with the pleiotropy test is its sensitivity to
start and end points for the two-dimensional scan. If the scan interval
is too small, the test can return misleading results. If *Tmem189* were
to have a very broad univariate peak, such that the peak extended beyond
the scan interval, then we might observe misleading results in tests
involving *Tmem189*.

``` r
# read oddball points from interactive plot
tr <- c("Tmem189",
        "Myo15b",
        "Sephs2",
        "Bcmo1",
        "Tmprss4",
        "Cacnb3",
        "Gm12929",
        "Eya2"
        )
tr_tib <- hot_annot %>%
  filter(symbol %in% tr) %>%
  arrange(annot.id)
K <- readRDS(here("analysis", "data", "derived_data", "kinship.rds"))
map <- readRDS(here("analysis", "data", "derived_data", "map.rds"))
covar <- readRDS(here("analysis", "data", "derived_data", "covar.rds"))
aprobs <- readRDS(here("analysis", "data", "derived_data", "genoprobs.rds"))
expr_sm <- expr[ , colnames(expr) %in% tr_tib$annot.id, drop = FALSE]
s1 <- qtl2::scan1(genoprobs = aprobs,
            pheno = expr_sm,
            kinship = K,
            addcovar = covar, 
            cores = 0
)
# zoom in on Chr2 after 150 Mb
# get length of map$`1`
purrr::map(.x = 1:8, .f = function(x) qtl2::plot_scan1(s1, map, lodcolumn = x, chr = 2, main = paste0(tr_tib[x, ])))
qtl2::find_peaks(scan1_output = s1, map = map, threshold = 7) %>% 
  as_tibble() %>%
  filter(chr == 2) %>%
  left_join(hot_annot, by = c("lodcolumn" = "annot.id"))
```

``` r
hot_chr <- 5
hot_mid <- 146
keller_mediator <- "Pdx1"
```

# Chr 5

``` r
hot_lower <- hot_mid - 2
hot_upper <- hot_mid + 2
fp <- paste0("Chr", hot_chr, "-")
knitr::opts_chunk$set(fig.path = here("analysis", "figures", fp))
```

``` r
pleio_mixed_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio.rds")))
pleio_local_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio_local.rds")))
pleio_nonlocal_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio_nonlocal.rds")))
med_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_mediations.rds")))
```

``` r
local_annot <- lod_peaks %>%
  filter(chrom == hot_chr, pos <= hot_upper, pos >= hot_lower) %>%
  left_join(y = annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= ult) %>%
  mutate(local_tx = chr == hot_chr) %>%
  filter(local_tx)
```

``` r
nonlocal_annot <- lod_peaks %>%
  filter(chrom == hot_chr, pos <= hot_upper, pos >= hot_lower) %>%
  inner_join(y = annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= 7.18) %>%
  mutate(local_tx = chr == hot_chr) %>%
  filter(!local_tx) %>%
  filter(!is.na(hotspot))
```

``` r
local_annot2 <- local_annot %>%
  filter(middle <= hot_upper, middle >= hot_lower) %>% 
  arrange(desc(lod))
```

``` r
hot_local <- expr[, colnames(expr) %in% local_annot2$annot.id]
hot_nonlocal <- expr[, colnames(expr) %in% nonlocal_annot$annot.id]
nnonlocal <- nrow(nonlocal_annot)
nlocal <- nrow(local_annot2)
```

``` r
med_tib <- med_list %>%
  purrr::map(bind_rows) %>%
  bind_rows() %>%
  as_tibble() %>%
  rename(lod_without_med = lod_no_med, lod_with_med = lod_med) %>%
  mutate(lod_diff = lod_without_med - lod_with_med, lod_diff_proportion = lod_diff / lod_without_med) %>%
  left_join(local_annot2, by = c("local_gene_id" = "annot.id")) %>%
  left_join(nonlocal_annot, by = c("nonlocal_gene_id" = "annot.id")) %>%
  dplyr::select(symbol.x, symbol.y, lod_without_med, lod_with_med, lod_diff, lod_diff_proportion, local_gene_id, nonlocal_gene_id) %>%
  rename(local_symbol = symbol.x, nonlocal_symbol = symbol.y) %>%
  filter(!is.na(local_symbol))
```

``` r
hot_annot <- local_annot2 %>%
  bind_rows(nonlocal_annot) %>%
  dplyr::select(annot.id, symbol)
pleio_tib <- pleio_mixed_list %>%
  purrr::map(.f = function(x) rename(x, gene1_id = local_gene_id, gene2_id = nonlocal_gene_id)) %>%
  c(pleio_local_list, pleio_nonlocal_list) %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(hot_annot, by = c("gene1_id" = "annot.id")) %>%
  rename(gene1_symbol = symbol) %>%
  inner_join(hot_annot, by = c("gene2_id" = "annot.id")) %>%
  rename(gene2_symbol = symbol) %>%
  dplyr::select(gene1_symbol, gene2_symbol, pleiotropy_lod)
```

``` r
pm <- pleio_tib %>%
  tibble_to_matrix(symmetric = TRUE)
```

``` r
p <- pleio_mixed_list %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(med_tib, by = c("local_gene_id", "nonlocal_gene_id")) %>%
  dplyr::select(local_symbol, nonlocal_symbol, pleiotropy_lod, lod_diff, lod_diff_proportion, lod_without_med, lod_with_med, filename, nonlocal_gene_id, local_gene_id) %>%
  ggplot() + geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = local_symbol == keller_mediator, nonlocal = nonlocal_symbol, local = local_symbol), size = 0.05) + theme(legend.position = "none") + xlab("Pleiotropy LOD") + ylab("Mediation LOD difference")
```

    ## Warning: Ignoring unknown aesthetics: nonlocal, local

``` r
fnsvg <- paste0("Chr", hot_chr, "_scatter.svg")
fnpdf <- paste0("Chr", hot_chr, "_scatter.pdf")
ggsave(plot = p, filename = here("analysis", "figures", fnsvg))
```

    ## Saving 7 x 5 in image

``` r
ggsave(plot = p, filename = here("analysis", "figures", fnpdf))
```

    ## Saving 7 x 5 in image

``` r
plotly::ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr5-unnamed-chunk-38-1.png)<!-- -->

``` r
foo <- pleio_mixed_list %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(med_tib, by = c("local_gene_id", "nonlocal_gene_id")) %>%
  dplyr::select(local_symbol, nonlocal_symbol, pleiotropy_lod, lod_diff, lod_diff_proportion, lod_without_med, lod_with_med, filename, nonlocal_gene_id, local_gene_id)
p <- foo %>%
  ggplot() +
  geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = nonlocal_symbol), size = 0.1) +
  facet_wrap(. ~ local_symbol, nrow = ceiling(nlocal / 2), ncol = 2) +
  broman::karl_theme(strip.placement = "inside") +
  theme(legend.position = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
) +
  xlab("Pleiotropy LOD") +
  ylab("Mediation LOD difference") +
  geom_text(aes(x = max(foo$pleiotropy_lod), y = 1, label = local_symbol), inherit.aes=FALSE)

fnsvg <- paste0("Chr", hot_chr, "_scatter-panel.svg")
fnpdf <- paste0("Chr", hot_chr, "_scatter-panel.pdf")
ggsave(p, filename = here("analysis", "figures", fnsvg), height = ceiling(nlocal / 2))
```

    ## Saving 7 x 11 in image

``` r
ggsave(p, filename = here("analysis", "figures", fnpdf), height = ceiling(nlocal / 2))
```

    ## Saving 7 x 11 in image

``` r
plotly::ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr5-unnamed-chunk-39-1.png)<!-- -->

``` r
# put in multiple figures, if there are many local traits
n_figs <- 1 + nlocal %/% 18 # 18 is the number of panels per page in pdf
n_last <- nlocal %% 18 # number of panels on last figure
max_pleio_lod <- max(foo$pleiotropy_lod)
max_med_lod_difference <- max(foo$lod_diff)
min_med_lod_difference <- min(foo$lod_diff)

for (i in 1:n_figs){
  n_panels <- 18 * (i < n_figs) + n_last * (i == n_figs)
  foo %>%
    filter(local_symbol %in% local_annot2$symbol[((i - 1) * 18 + 1): (18 * i)]) %>%
    ggplot() +
      geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = nonlocal_symbol), size = 0.1) +
      facet_wrap(. ~ local_symbol, nrow = ceiling(n_panels / 2), ncol = 2) +
      broman::karl_theme() +
      theme(legend.position = "none") +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
      xlab("Pleiotropy LOD") + xlim(c(0, max_pleio_lod)) + ylim(c(min_med_lod_difference, max_med_lod_difference)) +
      ylab("Mediation LOD difference") +
      geom_text(size = 2, aes(x = max_pleio_lod / 2, y = 2 * max_med_lod_difference / 3, label = local_symbol), inherit.aes=FALSE)
  fnsvg <- paste0("Chr", hot_chr, "_scatter-panel", i, ".svg")
  fnpdf <- paste0("Chr", hot_chr, "_scatter-panel", i, ".pdf")
  ggsave(filename = here("analysis", "figures", fnsvg), height = ceiling(n_panels / 2))
  ggsave(filename = here("analysis", "figures", fnpdf), height = ceiling(n_panels / 2))
}
```

``` r
scols1 <- 1L + (rownames(pm) %in% local_annot2$symbol)
scols2 <- rownames(pm) %in% (med_tib %>% filter(local_symbol == keller_mediator, lod_diff > 1.5))$nonlocal_symbol
```

``` r
# order the nonlocal traits first
## indicator for being a nonlocal trait
indic_nonlocal <- rownames(pm) %in% nonlocal_annot$symbol
# define ord_nonlocal
ord_nonlocal <- pm[indic_nonlocal, indic_nonlocal] %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(pm[indic_nonlocal, indic_nonlocal])) %>%
  order.dendrogram()
rn_nonlocal <- rownames(pm[indic_nonlocal, indic_nonlocal][ord_nonlocal, ord_nonlocal])
# define ord_local
ord_local<- pm[!indic_nonlocal, !indic_nonlocal] %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(pm[!indic_nonlocal, !indic_nonlocal])) %>%
  order.dendrogram()
rn_local <- rownames(pm[!indic_nonlocal, !indic_nonlocal][ord_local, ord_local])
rn <- c(rn_nonlocal, rn_local)
```

``` r
scols_extra <- purrr::map(.x = local_annot2$symbol,
                          .f = function(x) {rownames(pm) %in% (med_tib %>% filter(local_symbol == x, lod_diff > 1.5))$nonlocal_symbol}
) %>%
  bind_cols() %>%
  (function(x){colnames(x) <- local_annot2$symbol; x}) %>%
  mutate(local = scols1)
```

``` r
n_med_tr <- 4 # arbitrary threshold for number of mediated traits per local
# only those local traits with n_med_tr or more mediated nonlocal traits 
# appear as annotation columns in the heatmap
n_med_per_local <- med_tib %>%
  group_by(local_symbol) %>%
  filter(lod_diff > 1.5) %>%
  count() %>%
  filter(n >= n_med_tr)
n_med_per_local %>%
  arrange(desc(n)) %>%
  knitr::kable()
```

| local\_symbol |   n |
| :------------ | --: |
| Rpl21         | 105 |
| Slc46a3       |  96 |
| Polr1d        |  85 |
| Pdx1          |  74 |
| Tmem130       |  74 |
| Pan3          |  44 |
| Gm27033       |  37 |
| Gtf3a         |  33 |
| Cpsf4         |  22 |
| Arpc1a        |   9 |
| 2210019I11Rik |   7 |
| Rnf6          |   7 |
| Zscan25       |   6 |
| Nptx2         |   5 |

``` r
scols_pre <- scols_extra %>%
  dplyr::select(n_med_per_local$local_symbol)
# get ordering for columns of scols_pre
col_order <- scols_pre %>%
  t() %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(t(scols_pre))) %>%
  order.dendrogram()
scols <- scols_pre[ , col_order] %>%
  mutate(local = scols1)
fnhtml <- here("analysis", "figures", paste0("Chr", hot_chr, "hm.html"))
```

``` r
ord <- match(rn, rownames(pm))

(p <- pm[ord, ord] %>%
  apply(FUN = trunc2, MARGIN = 2, threshold = 5) %>%
  heatmaply(Rowv = FALSE,
            Colv = FALSE,
            cexRow = 0.2,
            cexCol = 0.2,
            symm = TRUE,
            row_side_colors = scols[ord, ],
            col_side_colors = tibble(local = scols1[ord]),
            limits = c(0, 5)
            ) %>%
  hide_guides()
)
```

    ## Warning in heatmaply.heatmapr(hm, colors = colors, limits = limits,
    ## scale_fill_gradient_fun = scale_fill_gradient_fun, : The hover text for
    ## col_side_colors is currently not implemented (due to an issue in plotly).
    ## We hope this would get resolved in future releases.

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr5-unnamed-chunk-45-1.png)<!-- -->

``` r
fnpdf <- here("analysis", "figures", paste0("Chr", hot_chr, "hm.pdf"))
htmlwidgets::saveWidget(widget = p, 
                        file = fnhtml
                        )
```

``` r
webshot::webshot(url = fnhtml,
                 file = fnpdf,
                 cliprect = "viewport"
                 )
```

``` r
bar <- foo %>%
  filter(local_symbol == keller_mediator)

p <- tibble(nonlocal_symbol = rownames(pm[indic_nonlocal, indic_nonlocal])[ord_nonlocal], index = 1:length(nonlocal_symbol)) %>%
  inner_join(bar, by = "nonlocal_symbol") %>% 
  left_join(nonlocal_annot, by = c("nonlocal_symbol" = "symbol")) %>%
  ggplot() + geom_point(aes(y = pos, x = index, color = lod_diff > 1.5, nonlocal = nonlocal_symbol, lod_diff = lod_diff)) 
```

    ## Warning: Ignoring unknown aesthetics: nonlocal, lod_diff

``` r
ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr5-unnamed-chunk-47-1.png)<!-- -->

## Chr 5 results

``` r
hot_chr <- 7
hot_mid <- 46
keller_mediator <- "Fam83e"
```

# Chr 7

``` r
hot_lower <- hot_mid - 2
hot_upper <- hot_mid + 2
fp <- paste0("Chr", hot_chr, "-")
knitr::opts_chunk$set(fig.path = here("analysis", "figures", fp))
```

``` r
pleio_mixed_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio.rds")))
pleio_local_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio_local.rds")))
pleio_nonlocal_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_pleio_nonlocal.rds")))
med_list <- readRDS(here("analysis", "data", "derived_data", paste0("Chr", hot_chr, "_mediations.rds")))
```

``` r
local_annot <- lod_peaks %>%
  filter(chrom == hot_chr, pos <= hot_upper, pos >= hot_lower) %>%
  left_join(y = annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= ult) %>%
  mutate(local_tx = chr == hot_chr) %>%
  filter(local_tx)
```

``` r
nonlocal_annot <- lod_peaks %>%
  filter(chrom == hot_chr, pos <= hot_upper, pos >= hot_lower) %>%
  inner_join(y = annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= 7.18) %>%
  mutate(local_tx = chr == hot_chr) %>%
  filter(!local_tx) %>%
  filter(!is.na(hotspot))
```

``` r
local_annot2 <- local_annot %>%
  filter(middle <= hot_upper, middle >= hot_lower) %>% 
  arrange(desc(lod))
```

``` r
hot_local <- expr[, colnames(expr) %in% local_annot2$annot.id]
hot_nonlocal <- expr[, colnames(expr) %in% nonlocal_annot$annot.id]
nnonlocal <- nrow(nonlocal_annot)
nlocal <- nrow(local_annot2)
```

``` r
med_tib <- med_list %>%
  purrr::map(bind_rows) %>%
  bind_rows() %>%
  as_tibble() %>%
  rename(lod_without_med = lod_no_med, lod_with_med = lod_med) %>%
  mutate(lod_diff = lod_without_med - lod_with_med, lod_diff_proportion = lod_diff / lod_without_med) %>%
  left_join(local_annot2, by = c("local_gene_id" = "annot.id")) %>%
  left_join(nonlocal_annot, by = c("nonlocal_gene_id" = "annot.id")) %>%
  dplyr::select(symbol.x, symbol.y, lod_without_med, lod_with_med, lod_diff, lod_diff_proportion, local_gene_id, nonlocal_gene_id) %>%
  rename(local_symbol = symbol.x, nonlocal_symbol = symbol.y) %>%
  filter(!is.na(local_symbol))
```

``` r
hot_annot <- local_annot2 %>%
  bind_rows(nonlocal_annot) %>%
  dplyr::select(annot.id, symbol)
pleio_tib <- pleio_mixed_list %>%
  purrr::map(.f = function(x) rename(x, gene1_id = local_gene_id, gene2_id = nonlocal_gene_id)) %>%
  c(pleio_local_list, pleio_nonlocal_list) %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(hot_annot, by = c("gene1_id" = "annot.id")) %>%
  rename(gene1_symbol = symbol) %>%
  inner_join(hot_annot, by = c("gene2_id" = "annot.id")) %>%
  rename(gene2_symbol = symbol) %>%
  dplyr::select(gene1_symbol, gene2_symbol, pleiotropy_lod)
```

``` r
pm <- pleio_tib %>%
  tibble_to_matrix(symmetric = TRUE)
```

``` r
p <- pleio_mixed_list %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(med_tib, by = c("local_gene_id", "nonlocal_gene_id")) %>%
  dplyr::select(local_symbol, nonlocal_symbol, pleiotropy_lod, lod_diff, lod_diff_proportion, lod_without_med, lod_with_med, filename, nonlocal_gene_id, local_gene_id) %>%
  ggplot() + geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = local_symbol == keller_mediator, nonlocal = nonlocal_symbol, local = local_symbol), size = 0.05) + theme(legend.position = "none") + xlab("Pleiotropy LOD") + ylab("Mediation LOD difference")
```

    ## Warning: Ignoring unknown aesthetics: nonlocal, local

``` r
fnsvg <- paste0("Chr", hot_chr, "_scatter.svg")
fnpdf <- paste0("Chr", hot_chr, "_scatter.pdf")
ggsave(plot = p, filename = here("analysis", "figures", fnsvg))
```

    ## Saving 7 x 5 in image

``` r
ggsave(plot = p, filename = here("analysis", "figures", fnpdf))
```

    ## Saving 7 x 5 in image

``` r
plotly::ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr7-unnamed-chunk-57-1.png)<!-- -->

``` r
foo <- pleio_mixed_list %>%
  bind_rows() %>%
  rename(pleiotropy_lod = lod) %>%
  inner_join(med_tib, by = c("local_gene_id", "nonlocal_gene_id")) %>%
  dplyr::select(local_symbol, nonlocal_symbol, pleiotropy_lod, lod_diff, lod_diff_proportion, lod_without_med, lod_with_med, filename, nonlocal_gene_id, local_gene_id)
p <- foo %>%
  ggplot() +
  geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = nonlocal_symbol), size = 0.1) +
  facet_wrap(. ~ local_symbol, nrow = ceiling(nlocal / 2), ncol = 2) +
  broman::karl_theme(strip.placement = "inside") +
  theme(legend.position = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
) +
  xlab("Pleiotropy LOD") +
  ylab("Mediation LOD difference") +
  geom_text(aes(x = max(foo$pleiotropy_lod), y = 1, label = local_symbol), inherit.aes=FALSE)

fnsvg <- paste0("Chr", hot_chr, "_scatter-panel.svg")
fnpdf <- paste0("Chr", hot_chr, "_scatter-panel.pdf")
ggsave(p, filename = here("analysis", "figures", fnsvg), height = ceiling(nlocal / 2))
```

    ## Saving 7 x 22 in image

``` r
ggsave(p, filename = here("analysis", "figures", fnpdf), height = ceiling(nlocal / 2))
```

    ## Saving 7 x 22 in image

``` r
plotly::ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr7-unnamed-chunk-58-1.png)<!-- -->

``` r
# put in multiple figures, if there are many local traits
n_figs <- 1 + nlocal %/% 18 # 18 is the number of panels per page in pdf
n_last <- nlocal %% 18 # number of panels on last figure
max_pleio_lod <- max(foo$pleiotropy_lod)
max_med_lod_difference <- max(foo$lod_diff)
min_med_lod_difference <- min(foo$lod_diff)

for (i in 1:n_figs){
  n_panels <- 18 * (i < n_figs) + n_last * (i == n_figs)
  foo %>%
    filter(local_symbol %in% local_annot2$symbol[((i - 1) * 18 + 1): (18 * i)]) %>%
    ggplot() +
      geom_point(mapping = aes(x = pleiotropy_lod, y = lod_diff, colour = nonlocal_symbol), size = 0.1) +
      facet_wrap(. ~ local_symbol, nrow = ceiling(n_panels / 2), ncol = 2) +
      broman::karl_theme() +
      theme(legend.position = "none") +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
      xlab("Pleiotropy LOD") + xlim(c(0, max_pleio_lod)) + ylim(c(min_med_lod_difference, max_med_lod_difference)) +
      ylab("Mediation LOD difference") +
      geom_text(size = 2, aes(x = max_pleio_lod / 2, y = 2 * max_med_lod_difference / 3, label = local_symbol), inherit.aes=FALSE)
  fnsvg <- paste0("Chr", hot_chr, "_scatter-panel", i, ".svg")
  fnpdf <- paste0("Chr", hot_chr, "_scatter-panel", i, ".pdf")
  ggsave(filename = here("analysis", "figures", fnsvg), height = ceiling(n_panels / 2))
  ggsave(filename = here("analysis", "figures", fnpdf), height = ceiling(n_panels / 2))
}
```

``` r
scols1 <- 1L + (rownames(pm) %in% local_annot2$symbol)
scols2 <- rownames(pm) %in% (med_tib %>% filter(local_symbol == keller_mediator, lod_diff > 1.5))$nonlocal_symbol
```

``` r
# order the nonlocal traits first
## indicator for being a nonlocal trait
indic_nonlocal <- rownames(pm) %in% nonlocal_annot$symbol
# define ord_nonlocal
ord_nonlocal <- pm[indic_nonlocal, indic_nonlocal] %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(pm[indic_nonlocal, indic_nonlocal])) %>%
  order.dendrogram()
rn_nonlocal <- rownames(pm[indic_nonlocal, indic_nonlocal][ord_nonlocal, ord_nonlocal])
# define ord_local
ord_local<- pm[!indic_nonlocal, !indic_nonlocal] %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(pm[!indic_nonlocal, !indic_nonlocal])) %>%
  order.dendrogram()
rn_local <- rownames(pm[!indic_nonlocal, !indic_nonlocal][ord_local, ord_local])
rn <- c(rn_nonlocal, rn_local)
```

``` r
scols_extra <- purrr::map(.x = local_annot2$symbol,
                          .f = function(x) {rownames(pm) %in% (med_tib %>% filter(local_symbol == x, lod_diff > 1.5))$nonlocal_symbol}
) %>%
  bind_cols() %>%
  (function(x){colnames(x) <- local_annot2$symbol; x}) %>%
  mutate(local = scols1)
```

``` r
n_med_tr <- 4 # arbitrary threshold for number of mediated traits per local
# only those local traits with n_med_tr or more mediated nonlocal traits 
# appear as annotation columns in the heatmap
n_med_per_local <- med_tib %>%
  group_by(local_symbol) %>%
  filter(lod_diff > 1.5) %>%
  count() %>%
  filter(n >= n_med_tr)
n_med_per_local %>%
  arrange(desc(n)) %>%
  knitr::kable()
```

| local\_symbol |  n |
| :------------ | -: |
| Fam83e        | 96 |
| Dhdh          | 67 |
| Cpt1c         | 62 |
| Rras          | 53 |
| Fuz           | 52 |
| Irf3          | 48 |
| Emp3          | 41 |
| Myh14         | 41 |
| Fam71e1       | 40 |
| Aldh16a1      | 39 |
| Kcnj14        | 39 |
| Pnkp          | 38 |
| Gm15545       | 32 |
| Abcc8         | 31 |
| Tmem86a       | 30 |
| Klk1b4        | 29 |
| Ppfia3        | 28 |
| Uevld         | 22 |
| Nosip         | 21 |
| Mtag2         | 19 |
| Rps11         | 15 |
| Tph1          | 14 |
| Ush1c         | 12 |
| Pih1d1        | 11 |
| Plekha4       |  8 |
| Fut1          |  7 |
| Izumo1        |  7 |
| Ntn5          |  7 |
| Shank1        |  6 |
| Gm16022       |  4 |
| Saal1         |  4 |
| Tmem143       |  4 |

``` r
scols_pre <- scols_extra %>%
  dplyr::select(n_med_per_local$local_symbol)
# get ordering for columns of scols_pre
col_order <- scols_pre %>%
  t() %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  dendextend::seriate_dendrogram(dist(t(scols_pre))) %>%
  order.dendrogram()
scols <- scols_pre[ , col_order] %>%
  mutate(local = scols1)
fnhtml <- here("analysis", "figures", paste0("Chr", hot_chr, "hm.html"))
```

``` r
ord <- match(rn, rownames(pm))

(p <- pm[ord, ord] %>%
  apply(FUN = trunc2, MARGIN = 2, threshold = 5) %>%
  heatmaply(Rowv = FALSE,
            Colv = FALSE,
            cexRow = 0.2,
            cexCol = 0.2,
            symm = TRUE,
            row_side_colors = scols[ord, ],
            col_side_colors = tibble(local = scols1[ord]),
            limits = c(0, 5)
            ) %>%
  hide_guides()
)
```

    ## Warning in heatmaply.heatmapr(hm, colors = colors, limits = limits,
    ## scale_fill_gradient_fun = scale_fill_gradient_fun, : The hover text for
    ## col_side_colors is currently not implemented (due to an issue in plotly).
    ## We hope this would get resolved in future releases.

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr7-unnamed-chunk-64-1.png)<!-- -->

``` r
fnpdf <- here("analysis", "figures", paste0("Chr", hot_chr, "hm.pdf"))
htmlwidgets::saveWidget(widget = p, 
                        file = fnhtml
                        )
```

``` r
webshot::webshot(url = fnhtml,
                 file = fnpdf,
                 cliprect = "viewport"
                 )
```

``` r
bar <- foo %>%
  filter(local_symbol == keller_mediator)

p <- tibble(nonlocal_symbol = rownames(pm[indic_nonlocal, indic_nonlocal])[ord_nonlocal], index = 1:length(nonlocal_symbol)) %>%
  inner_join(bar, by = "nonlocal_symbol") %>% 
  left_join(nonlocal_annot, by = c("nonlocal_symbol" = "symbol")) %>%
  ggplot() + geom_point(aes(y = pos, x = index, color = lod_diff > 1.5, nonlocal = nonlocal_symbol, lod_diff = lod_diff)) 
```

    ## Warning: Ignoring unknown aesthetics: nonlocal, lod_diff

``` r
ggplotly(p)
```

![](/Users/frederickboehm/Box%20Sync/attie/qtl2hotspots/analysis/figures/Chr7-unnamed-chunk-66-1.png)<!-- -->

``` r
hot_chr <- 11
hot_mid <- 71
keller_mediator <- "Sat2"
```

``` r
hot_chr <- 13
hot_mid <- 112.5
keller_mediator <- "Il6st"
```

# Discussion

Pleiotropy tests complement mediation analyses in the following manner.
When considering a collection of potential mediators for a nonlocal
trait, mediation analyses clearly identify the mediator when there is,
in fact, a mediator among the candidates. However, when there is no
mediator, mediation analyses provide little information about the
genetics of the complex trait. Pleiotropy testing is particularly useful
in this setting, since a pleiotropy test always gives an answer about
the number of underlying QTL. A drawback of pleiotropy testing is that
the scientific question that it addresses - do two traits share a single
QTL - informs, but doesn’t fully resolve, the genetics of the nonlocal
complex traits under study.

# Colophon

This report was generated on `r Sys.time()` using the following
computational environment and dependencies:

``` r
# which R packages and versions?
devtools::session_info()
```

    ## ─ Session info ──────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 3.6.1 (2019-07-05)
    ##  os       macOS Mojave 10.14.6        
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2019-07-31                  
    ## 
    ## ─ Packages ──────────────────────────────────────────────────────────────
    ##  package      * version  date       lib
    ##  assertthat     0.2.1    2019-03-21 [1]
    ##  backports      1.1.4    2019-04-10 [1]
    ##  bitops         1.0-6    2013-08-17 [1]
    ##  bookdown       0.12     2019-07-11 [1]
    ##  broman       * 0.69-5   2019-04-11 [1]
    ##  broom          0.5.2    2019-04-07 [1]
    ##  callr          3.3.1    2019-07-18 [1]
    ##  caTools        1.17.1.2 2019-03-06 [1]
    ##  cellranger     1.1.0    2016-07-27 [1]
    ##  cli            1.1.0    2019-03-19 [1]
    ##  cluster        2.1.0    2019-06-19 [1]
    ##  codetools      0.2-16   2018-12-24 [1]
    ##  colorspace     1.4-1    2019-03-18 [1]
    ##  crayon         1.3.4    2017-09-16 [1]
    ##  crosstalk      1.0.0    2016-12-21 [1]
    ##  data.table     1.12.2   2019-04-07 [1]
    ##  dendextend     1.12.0   2019-05-11 [1]
    ##  desc           1.2.0    2018-05-01 [1]
    ##  devtools       2.1.0    2019-07-06 [1]
    ##  digest         0.6.20   2019-07-04 [1]
    ##  dplyr        * 0.8.3    2019-07-04 [1]
    ##  evaluate       0.14     2019-05-28 [1]
    ##  forcats      * 0.4.0    2019-02-17 [1]
    ##  foreach        1.4.7    2019-07-27 [1]
    ##  fs             1.3.1    2019-05-06 [1]
    ##  gclus          1.3.2    2019-01-07 [1]
    ##  gdata          2.18.0   2017-06-06 [1]
    ##  gdtools      * 0.1.9    2019-06-18 [1]
    ##  generics       0.0.2    2018-11-29 [1]
    ##  ggplot2      * 3.2.0    2019-06-16 [1]
    ##  git2r          0.26.1   2019-06-29 [1]
    ##  glue           1.3.1    2019-03-12 [1]
    ##  gplots         3.0.1.1  2019-01-27 [1]
    ##  gridExtra      2.3      2017-09-09 [1]
    ##  gtable         0.3.0    2019-03-25 [1]
    ##  gtools         3.8.1    2018-06-26 [1]
    ##  haven          2.1.1    2019-07-04 [1]
    ##  heatmaply    * 1.0.0    2019-07-26 [1]
    ##  here         * 0.1      2017-05-28 [1]
    ##  highr          0.8      2019-03-20 [1]
    ##  hms            0.5.0    2019-07-09 [1]
    ##  htmltools      0.3.6    2017-04-28 [1]
    ##  htmlwidgets    1.3      2018-09-30 [1]
    ##  httpuv         1.5.1    2019-04-05 [1]
    ##  httr           1.4.0    2018-12-11 [1]
    ##  iterators      1.0.12   2019-07-26 [1]
    ##  jsonlite       1.6      2018-12-07 [1]
    ##  KernSmooth     2.23-15  2015-06-29 [1]
    ##  knitr          1.23     2019-05-18 [1]
    ##  labeling       0.3      2014-08-23 [1]
    ##  later          0.8.0    2019-02-11 [1]
    ##  lattice        0.20-38  2018-11-04 [1]
    ##  lazyeval       0.2.2    2019-03-15 [1]
    ##  lubridate      1.7.4    2018-04-11 [1]
    ##  magrittr       1.5      2014-11-22 [1]
    ##  MASS           7.3-51.4 2019-03-31 [1]
    ##  memoise        1.1.0    2017-04-21 [1]
    ##  mime           0.7      2019-06-11 [1]
    ##  modelr         0.1.4    2019-02-18 [1]
    ##  munsell        0.5.0    2018-06-12 [1]
    ##  nlme           3.1-140  2019-05-12 [1]
    ##  pillar         1.4.2    2019-06-29 [1]
    ##  pkgbuild       1.0.3    2019-03-20 [1]
    ##  pkgconfig      2.0.2    2018-08-16 [1]
    ##  pkgload        1.0.2    2018-10-29 [1]
    ##  plotly       * 4.9.0    2019-04-10 [1]
    ##  plyr           1.8.4    2016-06-08 [1]
    ##  prettyunits    1.0.2    2015-07-13 [1]
    ##  processx       3.4.1    2019-07-18 [1]
    ##  promises       1.0.1    2018-04-13 [1]
    ##  ps             1.3.0    2018-12-21 [1]
    ##  purrr        * 0.3.2    2019-03-15 [1]
    ##  qtl2hotspots * 0.1.0    2019-07-31 [1]
    ##  R6             2.4.0    2019-02-14 [1]
    ##  RColorBrewer   1.1-2    2014-12-07 [1]
    ##  Rcpp           1.0.2    2019-07-25 [1]
    ##  readr        * 1.3.1    2018-12-21 [1]
    ##  readxl         1.3.1    2019-03-13 [1]
    ##  registry       0.5-1    2019-03-05 [1]
    ##  remotes        2.1.0    2019-06-24 [1]
    ##  reshape2       1.4.3    2017-12-11 [1]
    ##  rlang          0.4.0    2019-06-25 [1]
    ##  rmarkdown      1.14     2019-07-12 [1]
    ##  rprojroot      1.3-2    2018-01-03 [1]
    ##  rstudioapi     0.10     2019-03-19 [1]
    ##  rvest          0.3.4    2019-05-15 [1]
    ##  scales         1.0.0    2018-08-09 [1]
    ##  seriation      1.2-7    2019-06-08 [1]
    ##  sessioninfo    1.1.1    2018-11-05 [1]
    ##  shiny          1.3.2    2019-04-22 [1]
    ##  stringi        1.4.3    2019-03-12 [1]
    ##  stringr      * 1.4.0    2019-02-10 [1]
    ##  svglite        1.2.2    2019-05-17 [1]
    ##  testthat       2.2.1    2019-07-25 [1]
    ##  tibble       * 2.1.3    2019-06-06 [1]
    ##  tidyr        * 0.8.3    2019-03-01 [1]
    ##  tidyselect     0.2.5    2018-10-11 [1]
    ##  tidyverse    * 1.2.1    2017-11-14 [1]
    ##  TSP            1.1-7    2019-05-22 [1]
    ##  usethis        1.5.1    2019-07-04 [1]
    ##  vctrs          0.2.0    2019-07-05 [1]
    ##  viridis      * 0.5.1    2018-03-29 [1]
    ##  viridisLite  * 0.3.0    2018-02-01 [1]
    ##  webshot        0.5.1    2018-09-28 [1]
    ##  withr          2.1.2    2018-03-15 [1]
    ##  xfun           0.8      2019-06-25 [1]
    ##  xml2           1.2.0    2018-01-24 [1]
    ##  xtable         1.8-4    2019-04-21 [1]
    ##  yaml           2.2.0    2018-07-25 [1]
    ##  zeallot        0.1.0    2018-01-28 [1]
    ##  source                              
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  Github (talgalili/heatmaply@a352dff)
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  local                               
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.1)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ##  CRAN (R 3.6.0)                      
    ## 
    ## [1] /Library/Frameworks/R.framework/Versions/3.6/Resources/library

The current Git commit details are:

``` r
# what commit is this file at?
git2r::repository(here::here())
```

    ## Local:    master /Users/frederickboehm/Box Sync/attie/qtl2hotspots
    ## Remote:   master @ origin (https://github.com/fboehm/qtl2hotspots.git)
    ## Head:     [7a7dbb7] 2019-07-31: deleted makefile for pdf / Rnw

<div id="refs" class="references">

<div id="ref-gemma2">

Boehm, Frederick. 2018. *Gemma2: Zhou & Stephens (2014) Gemma
Multivariate Linear Mixed Model*. <https://github.com/fboehm/gemma2>.

</div>

<div id="ref-boehm2019thesis">

Boehm, Frederick J. 2019. “Testing Pleiotropy Vs. Separate Qtl in
Multiparental Populations.” PhD thesis, University of Wisconsin-Madison.

</div>

<div id="ref-boehm2019testing">

Boehm, Frederick J., Elissa J. Chesler, Brian S. Yandell, and Karl W.
Broman. 2019. “Testing Pleiotropy Vs. Separate Qtl in Multiparental
Populations.” *G3: Genes| Genomes| Genetics*. Genetics Soc America.

</div>

<div id="ref-broman2009guide">

Broman, Karl W., and Saunak Sen. 2009. *A Guide to Qtl Mapping with
R/Qtl*. Vol. 46. Springer.

</div>

<div id="ref-chick2016defining">

Chick, Joel M, Steven C Munger, Petr Simecek, Edward L Huttlin, Kwangbom
Choi, Daniel M Gatti, Narayanan Raghupathy, Karen L Svenson, Gary A
Churchill, and Steven P Gygi. 2016. “Defining the Consequences of
Genetic Variation on a Proteome-Wide Scale.” *Nature* 534 (7608). Nature
Publishing Group: 500.

</div>

<div id="ref-churchill2004collaborative">

Churchill, Gary A, David C Airey, Hooman Allayee, Joe M Angel, Alan D
Attie, Jackson Beatty, William D Beavis, et al. 2004. “The Collaborative
Cross, a Community Resource for the Genetic Analysis of Complex Traits.”
*Nature Genetics* 36 (11). Nature Publishing Group: 1133–7.

</div>

<div id="ref-dendextend">

Galili, Tal. 2015. “Dendextend: An R Package for Visualizing, Adjusting,
and Comparing Trees of Hierarchical Clustering.” *Bioinformatics*.
<https://doi.org/10.1093/bioinformatics/btv428>.

</div>

<div id="ref-hahsler2008getting">

Hahsler, Michael, Kurt Hornik, and Christian Buchta. 2008. “Getting
Things in Order: An Introduction to the R Package Seriation.” *Journal
of Statistical Software* 25 (3). American Statistical Association: 1–34.

</div>

<div id="ref-huang2012multiparent">

Huang, Bevan E, Andrew W George, Kerrie L Forrest, Andrzej Kilian,
Matthew J Hayden, Matthew K Morell, and Colin R Cavanagh. 2012. “A
Multiparent Advanced Generation Inter-Cross Population for Genetic
Analysis in Wheat.” *Plant Biotechnology Journal* 10 (7). Wiley Online
Library: 826–39.

</div>

<div id="ref-huang2011analysis">

Huang, Xueqing, Maria-Joao Paulo, Martin Boer, Sigi Effgen, Paul Keizer,
Maarten Koornneef, and Fred A van Eeuwijk. 2011. “Analysis of Natural
Allelic Variation in Arabidopsis Using a Multiparent Recombinant Inbred
Line Population.” *Proceedings of the National Academy of Sciences* 108
(11). National Acad Sciences: 4488–93.

</div>

<div id="ref-jansen2007quantitative">

Jansen, Ritsert C. 2007. “Quantitative Trait Loci in Inbred Lines.” In
*Handbook of Statistical Genetics*, edited by David J Balding, Martin
Bishop, and Chris Cannings, 587–622. John Wiley & Sons.

</div>

<div id="ref-jiang1995multiple">

Jiang, Changjian, and Zhao-Bang Zeng. 1995. “Multiple Trait Analysis of
Genetic Mapping for Quantitative Trait Loci.” *Genetics* 140 (3).
Genetics Soc America: 1111–27.

</div>

<div id="ref-keller2018genetic">

Keller, Mark P, Daniel M Gatti, Kathryn L Schueler, Mary E Rabaglia,
Donnie S Stapleton, Petr Simecek, Matthew Vincent, et al. 2018. “Genetic
Drivers of Pancreatic Islet Function.” *Genetics*. Genetics Soc America,
genetics–300864.

</div>

<div id="ref-de2017back">

Koning, Dirk-Jan de, and Lauren M McIntyre. 2017. “Back to the Future:
Multiparent Populations Provide the Key to Unlocking the Genetic Basis
of Complex Traits.” *Genetics* 206. Genetics Soc America: 527–29.

</div>

<div id="ref-kover2009multiparent">

Kover, Paula X, William Valdar, Joseph Trakalo, Nora Scarcelli, Ian M
Ehrenreich, Michael D Purugganan, Caroline Durrant, and Richard Mott.
2009. “A Multiparent Advanced Generation Inter-Cross to Fine-Map
Quantitative Traits in Arabidopsis Thaliana.” *PLoS Genetics* 5 (7).
Public Library of Science: e1000551.

</div>

<div id="ref-lance1967general">

Lance, Godfrey N, and William Thomas Williams. 1967. “A General Theory
of Classificatory Sorting Strategies: 1. Hierarchical Systems.” *The
Computer Journal* 9 (4). The British Computer Society: 373–80.

</div>

<div id="ref-lander1989mapping">

Lander, Eric S, and David Botstein. 1989. “Mapping Mendelian Factors
Underlying Quantitative Traits Using Rflp Linkage Maps.” *Genetics* 121
(1). Genetics Soc America: 185–99.

</div>

<div id="ref-r">

R Core Team. 2018. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

<div id="ref-sax1923association">

Sax, Karl. 1923. “The Association of Size Differences with Seed-Coat
Pattern and Pigmentation in Phaseolus Vulgaris.” *Genetics* 8 (6).
Genetics Society of America: 552.

</div>

<div id="ref-shivakumar2018soybean">

Shivakumar, M, Gireesh C Giriraj Kumawat, SV Ramesh, and SM Husain.
2018. “Soybean Magic Population: A Novel Resource for Genetics and Plant
Breeding.” Current Science.

</div>

<div id="ref-soller1976power">

Soller, M, T Brody, and A Genizi. 1976. “On the Power of Experimental
Designs for the Detection of Linkage Between Marker Loci and
Quantitative Loci in Crosses Between Inbred Lines.” *Theoretical and
Applied Genetics* 47 (1). Springer: 35–39.

</div>

<div id="ref-stanley2017genetic">

Stanley, Patrick D, Enoch Ngoma, Siri ODay, and Elizabeth G King. 2017.
“Genetic Dissection of Nutrition-Induced Plasticity in
Insulin/Insulin-Like Growth Factor Signaling and Median Life Span in a
Drosophila Multiparent Population.” *Genetics* 206 (2). Genetics Soc
America: 587–602.

</div>

<div id="ref-svenson2012high">

Svenson, Karen L, Daniel M Gatti, William Valdar, Catherine E Welsh,
Riyan Cheng, Elissa J Chesler, Abraham A Palmer, Leonard McMillan, and
Gary A Churchill. 2012. “High-Resolution Genetic Mapping Using the Mouse
Diversity Outbred Population.” *Genetics* 190 (2). Genetics Soc America:
437–47.

</div>

<div id="ref-tisne2017identification">

Tisne, Sebastien, Virginie Pomies, Virginie Riou, Indra Syahputra,
Benoit Cochard, and Marie Denis. 2017. “Identification of Ganoderma
Disease Resistance Loci Using Natural Field Infection of an Oil Palm
Multiparental Population.” *G3: Genes, Genomes, Genetics* 7 (6). G3:
Genes, Genomes, Genetics: 1683–92.

</div>

<div id="ref-ggplot2">

Wickham, Hadley. 2016. *Ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

</div>

<div id="ref-zhou2014efficient">

Zhou, Xiang, and Matthew Stephens. 2014. “Efficient Multivariate Linear
Mixed Model Algorithms for Genome-Wide Association Studies.” *Nature
Methods* 11 (4). Nature Research: 407–9.

</div>

</div>
