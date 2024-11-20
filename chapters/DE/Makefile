BACKEND = bibtex
output = DE
ENGINE = xelatex
$(output).pdf: $(output).Rnw
	Rscript -e "Sweave(commandArgs(TRUE))" $< 
	perl -0777 -i -pe "s/(?<=backend=)\s*\w+/$(BACKEND)/" $(output).tex
	$(ENGINE) $(output).tex
	$(BACKEND) $(output) && $(ENGINE) $(output).tex
	rm -f *-blx* *.bbl *.aux *.toc *.lof *.blg *.bcf *.log *.xml *.out *.lot Makeoutput.tex *.bcf *.nav *.snm

clean:
	rm -f $(TARGET) $(output).pdf *.tex .Rhistory
	rm -f -r *pdf *tex *log *aux *out *blx* *lot *toc *xml *blg *bbl *lof *bcf *nav *snm

