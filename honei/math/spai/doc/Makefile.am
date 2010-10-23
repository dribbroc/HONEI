
DOCNAME = spaidoc

EXTRA_DIST = $(DOCNAME).tex $(DOCNAME).ps $(DOCNAME).pdf

docdir   = $(prefix)/share/doc/spai
doc_DATA = $(DOCNAME).ps $(DOCNAME).pdf


%.ps: %.tex
	-latex --interaction nonstopmode $*
	latex --interaction nonstopmode $*
	dvips $*

%.pdf: %.ps
	-pdflatex --interaction nonstopmode $*
	pdflatex --interaction nonstopmode $*


%.html: %.tex clean
	perl -pi.HV -e 's/\\documentclass(\[.*\])?{.*}/\\documentclass{article}/g' $(DOCNAME).tex
	hevea $<
	hevea $<
	mv $(DOCNAME).tex.HV $(DOCNAME).tex


WWW_DIR = /mnt/Enterprise/WebServer/Computational/Documents/software/spai

www: $(DOCNAME).tex $(DOCNAME).html $(DOCNAME).ps $(DOCNAME).pdf
	cp $^ $(WWW_DIR)
	for file in $^; do  chmod a+r $(WWW_DIR)/$$file;  done


CLEANFILES = $(DOCNAME).haux $(DOCNAME).htoc				\
	$(DOCNAME).dvi $(DOCNAME).log $(DOCNAME).aux $(DOCNAME).toc	\
	$(DOCNAME).out

MAINTAINERCLEANFILES =							\
	$(DOCNAME).ps $(DOCNAME).pdf $(DOCNAME)*.html
