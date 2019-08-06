
PAPER_DIR=.
FIGURE_DIR=../../figures

# Rules

all: $(PAPER_DIR)/hotspots-paper.html

$(PAPER_DIR)/hotspots-paper.html: $(PAPER_DIR)/hotspots-paper.Rmd $(PAPER_DIR)/methods.Rmd $(PAPER_DIR)/results.Rmd $(PAPER_DIR)/intro.Rmd $(PAPER_DIR)/discussion.Rmd $(PAPER_DIR)/hotspot-one-chromosome.Rmd $(PAPER_DIR)/abstract.Rmd
	R -e "devtools::install(dep = TRUE)"; cd $(PAPER_DIR); R -e "rmarkdown::render('$(notdir $<)')"




.PHONY: all

clean:
	rm -f $(PAPER_DIR)/hotspots-paper.html
	rm -f $(PAPER_DIR)/hotspots-paper.md
	rm -f $(PAPER_DIR)/hotspots-paper.tex

cleanfigures:
	rm -f $(FIGURE_DIR)/*.pdf $(FIGURE_DIR)/*.html $(FIGURE_DIR)/*.svg; rm -rf $(FIGURE_DIR)/*_files

cleanall: clean cleanfigures

