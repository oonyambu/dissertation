EXT = .Rnw
BACKEND = bibtex
ENGINE = xelatex
chapter2 = chapters/DE/DE
chapter3 = chapters/EGO/EGO
chapter4 = chapters/RSO/RSO

include = $(wildcard include/*)

package = $(wildcard package/*)
output = $(basename $(wildcard *.Rnw))

TARGET = $(chapter2)_processed$(EXT) \
         $(chapter3)_processed$(EXT) \
         $(chapter4)_processed$(EXT)


pat1 = "s/^(?=.*(?=\n.title))|^(?=(?(?!\\\\title).)*(?=\n.section\\{Introduction\\}))|^(?=.end\\{document\\})|^(?=.printbibliography)/%/gms"
pat2 = "s/^(\\\\)title[{]\\s*(\\\\bf)?\\s*/\1chapter{/s"
pat3 = "s/(>>=)/\1\n dr <- getwd()\n setwd('"
pat4 = "s/(.*end\\{document\\})/\n<<fig=false,echo=false>>=\nsetwd(dr)\n@\n\n\1/s"
pat5 = "s/(includegraphics.*\\{)(pdfs\/)/"
pat6 = "s/^.section[*]?{Appendix: */\\\\appendix{/s"


	
all: $(output).pdf
	open $<

$(chapter2)_processed$(EXT): $(chapter2)$(EXT)
	perl -0777 -pe $(pat1) $< > $@
	perl -i -pe $(pat2) $@
	perl -0777 -i -pe $(pat3)$(subst /,\\/,$(dir $<))"')/s" $@
	perl -i -pe $(pat4) $@
	perl -i -pe $(pat5)"\1$(subst /,\\/,$(dir $<))\2/s" $@
	# perl -i -pe $(pat6) $@

$(chapter3)_processed$(EXT): $(chapter3)$(EXT)
	perl -0777 -pe $(pat1) $< > $@
	perl -i -pe $(pat2) $@
	perl -0777 -i -pe $(pat3)$(subst /,\\/,$(dir $<))"')/s" $@
	perl -i -pe $(pat4) $@
	perl -i -pe $(pat5)"\1$(subst /,\\/,$(dir $<))\2/s" $@
	# perl -i -pe $(pat6) $@

$(chapter4)_processed$(EXT): $(chapter4)$(EXT)
	perl -0777 -pe $(pat1) $< > $@
	perl -i -pe $(pat2) $@
	perl -0777 -i -pe $(pat3)$(subst /,\\/,$(dir $<))"')/s" $@
	perl -i -pe $(pat4) $@
	perl -i -pe $(pat5)"\1$(subst /,\\/,$(dir $<))\2/s" $@
	# perl -i -pe $(pat6) $@

$(output).tex: $(output).Rnw $(TARGET) $(include)
	Rscript -e "Sweave(commandArgs(TRUE))" $<
	

$(output).pdf: $(output).tex 
	perl -0777 -i -pe "s/(?<=backend=)\s*\w+/$(BACKEND)/" include/preamble.tex
	$(ENGINE) $(output).tex
	$(BACKEND) $(output) && $(ENGINE) $(output).tex
	rm -f *-blx* *.bbl *.aux *.toc *.lof *.blg *.bcf *.log *.xml *.out *.lot Makefile.tex
	perl -0777 -i -pe "s/(?<=backend=)\s*\w+/bibtex/" include/preamble.tex

clean:
	rm -f $(TARGET) $(output).pdf .Rhistory $(TARGET:_processed.Rnw=.tex)
	rm -f -r *pdf *tex *log *aux *out *blx* *lot *toc *xml *blg *bbl *lof 