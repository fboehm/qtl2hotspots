
PAPER_DIR=analysis/paper


# Rules

all: $(PAPER_DIR)/hotspots-paper.html

$(PAPER_DIR)/hotspots-paper.html: $(PAPER_DIR)/hotspots-paper.Rmd
	cd $(PAPER_DIR); R -e "rmarkdown::render('$(notdir $<)', output_format = 'all')"




.PHONY: all


