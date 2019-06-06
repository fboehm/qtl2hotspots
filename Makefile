

all: analysis/paper/paper.pdf analysis/paper/paper.html
.PHONY: all


analysis/paper/paper.pdf: analysis/paper/paper.Rmd analysis/paper/references.bib
  Rscript -e "rmarkdown::render('$<', output_format = 'all')"

analysis/paper/paper.html: analysis/paper/paper.Rmd analysis/paper/references.bib
  Rscript -e "rmarkdown::render('$<', output_format = 'all')"
