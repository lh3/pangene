.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:pangene.pdf

pangene.pdf:pangene.tex pangene.bib
		pdflatex pangene; bibtex pangene; pdflatex pangene; pdflatex pangene;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ pangene.pdf
