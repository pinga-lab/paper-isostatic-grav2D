help:
	@echo "Available makes:"
	@echo "manuscript clean"

all:
	make manuscript
	make figures
	make clean

# Compile the manuscript
manuscript:
	pdflatex manuscript.tex
	bibtex   bib-file
	pdflatex manuscript.tex
	pdflatex manuscript.tex

# Convert and rename figures
figures:
	

# Remove unused files
cleanall:
	rm *.aux *.bbl *.blg *.log *.pdf *.lof *.lot *.synctex.gz

