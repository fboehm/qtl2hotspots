
PAPER_DIR=analysis/paper

all: $(PAPER_DIR)/paper.pdf


$(PAPER_DIR)/paper.pdf: $(PAPER_DIR)/paper.Rmd
	R -e "devtools::install(dep = TRUE)"; Rscript -e "rmarkdown::render('$<', output_format = 'all')"

