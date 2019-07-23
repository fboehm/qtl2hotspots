
PAPER_DIR=analysis/paper
FIGURE_DIR=analysis/figures
# Tools
LATEXMK = latexmk
RM = rm -f

# Rules

all: $(PAPER_DIR)/paper.pdf


# use tinytex to find and install needed latex packages
$(PAPER_DIR)/paper.tex: $(PAPER_DIR)/paper.Rnw
	R -e "devtools::install(dep = TRUE)"; cd $(PAPER_DIR); Rscript -e "knitr::knit('$(notdir $<)')"

$(PAPER_DIR)/paper.pdf: $(PAPER_DIR)/paper.tex $(PAPER_DIR)/research.bib pdf_figures
	cd $(PAPER_DIR); R -e "install.packages('tinytex', repos = 'https://cloud.r-project.org')"; R -e "tinytex::latexmk('$(notdir $<)', bib_engine = 'biber')"

pdf_figures: $(PAPER_DIR)/figures.Rmd
	cd $(PAPER_DIR); R -e "rmarkdown::render('$(notdir $<)')"

mostlyclean:
	cd $(PAPER_DIR); $(LATEXMK) -silent -c
	$(RM) $(PAPER_DIR)/*.bbl

clean: mostlyclean
	cd $(PAPER_DIR); $(LATEXMK) -silent -C
	$(RM) $(PAPER_DIR)/*.run.xml $(PAPER_DIR)/*.synctex.gz
	$(RM) $(PAPER_DIR)/*.d
	$(RM) $(PAPER_DIR)/paper.tex
  $(RM) $(FIGURE_DIR)/*.*

.PHONY: all clean mostlyclean pdf_figures



# Include auto-generated dependencies
#-include *.d


# use latexdiff to get PDF comparing current to an older version
diff:
	latexdiff latexdiff/main_old.tex $(PAPER_DIR)/main.tex > latexdiff/diff.tex
	mv latexdiff/diff.tex $(PAPER_DIR)
	cd $(PAPER_DIR);$(LATEXMK) -pdfps -bibtex diff.tex
	mv $(PAPER_DIR)/diff* latexdiff/

