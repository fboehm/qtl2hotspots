
PAPER_DIR=analysis/paper


# Rules

all: $(PAPER_DIR)/hotspots-paper.html $(PAPER_DIR)/hotspots-paper.md

$(PAPER_DIR)/hotspots-paper.html: $(PAPER_DIR)/hotspots-paper.Rmd
	R -e "devtools::install(dep = TRUE)"; cd $(PAPER_DIR); R -e "rmarkdown::render('$(notdir $<)', output_format = 'all', output_options = list(pandoc_args = c('--filter', '/root/.cabal/bin/pandoc-crossref')))"

$(PAPER_DIR)/hotspots-paper.md: $(PAPER_DIR)/hotspots-paper.Rmd
	R -e "devtools::install(dep = TRUE)"; cd $(PAPER_DIR); R -e "rmarkdown::render('$(notdir $<)', output_format = 'all', output_options = list(pandoc_args = c('--filter', '/root/.cabal/bin/pandoc-crossref')))"



.PHONY: all

clean:
	rm -f $(PAPER_DIR)/hotspots-paper.html
	rm -f $(PAPER_DIR)/hotspots-paper.md
