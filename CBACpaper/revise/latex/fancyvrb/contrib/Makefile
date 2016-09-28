
# `Makefile' for `pst-3dplot.pdf', hv, 2007/03/17

.SUFFIXES : .tex .ltx .dvi .ps .pdf .eps

PACKAGE = fvrb-ex

MAIN = $(PACKAGE)

TDS = ~/PSTricks/PSTricks-TDS

LATEX = latex

ARCHNAME = $(MAIN)-$(shell date +%y%m%d)

ARCHFILES = $(PACKAGE).sty $(PACKAGE).tex $(MAIN).tex README Changes Makefile

all : doc clean 
doc: $(MAIN).pdf

$(MAIN).pdf : $(MAIN).ps
	GS_OPTIONS=-dAutoRotatePages=/None ps2pdf $<

$(MAIN).ps : $(MAIN).dvi
	dvips $<

$(MAIN).dvi : $(MAIN).dtx
	$(LATEX) $<
	makeindex -s gind.ist $(basename $<)
	makeindex -s gglo.ist -o $(basename $<).gls $(basename $<).glo
#	bibtex $(basename $<)
	$(LATEX) $<
	$(LATEX) $<

clean : 
	$(RM) $(addprefix $(MAIN), .cb .cb? .log .aux .glg .glo .gls .ilg .idx .ind .tmp .toc .out .blg .bbl )
	$(RM) $(addprefix $(MAIN), .dvi .ps .xcp)
	$(RM) file.txt hello.f90

veryclean : clean
	$(RM) $(addprefix $(MAIN), .pdf .bbl .blg)

arch :
	zip $(ARCHNAME).zip $(ARCHFILES)

tds:
	cp -u Changes     $(TDS)/doc/generic/$(PACKAGE)/
	cp -u README      $(TDS)/doc/generic/$(PACKAGE)/
	cp -u $(MAIN).pdf $(TDS)/doc/generic/$(PACKAGE)/
#
	cp -u Changes        $(TDS)/tex/latex/$(PACKAGE)/
	cp -u $(PACKAGE).sty $(TDS)/tex/latex/$(PACKAGE)/
#
	cp -u Changes        $(TDS)/tex/generic/$(PACKAGE)/
	cp -u $(PACKAGE).tex $(TDS)/tex/generic/$(PACKAGE)/
#
	cp -u Changes        $(TDS)/dvips/$(PACKAGE)/
	cp -u $(PACKAGE).pro $(TDS)/dvips/$(PACKAGE)/
	cp -u $(PACKAGE).pro ~/Links/dvips-local/
#
	cp -u Changes     $(TDS)/source/$(PACKAGE)/
	cp -u $(MAIN).tex $(TDS)/source/$(PACKAGE)/
	cp -u $(MAIN).bib $(TDS)/source/$(PACKAGE)/
	cp -u Makefile    $(TDS)/source/$(PACKAGE)/


# EOF
