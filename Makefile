
PAPER_DIR=analysis/paper


# Rules

all: $(PAPER_DIR)/hotspots-paper.html

$(PAPER_DIR)/hotspots-paper.html: $(PAPER_DIR)/hotspots-paper.Rmd $(PAPER_DIR)/methods.Rmd $(PAPER_DIR)/results.Rmd $(PAPER_DIR)/intro.Rmd $(PAPER_DIR)/hotspot-one-chromosome.Rmd
	R -e "devtools::install(dep = TRUE)"; cd $(PAPER_DIR); R -e "rmarkdown::render('$(notdir $<)', output_format = 'all')"




.PHONY: all

clean:
	rm -f $(PAPER_DIR)/hotspots-paper.html
