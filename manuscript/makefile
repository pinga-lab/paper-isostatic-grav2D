help:
	@echo "Available makes:"
	@echo "manuscript - compile the latex files to build the manuscript"
	@echo "figs - rename, convert to tiff and set colorpace, dpi and size of figures"
	@echo "clean - remove unused files"

all:
	make manuscript
	make figs
	make clean

# Compile the manuscript
manuscript:
	pdflatex manuscript.tex
	bibtex bib-file
	pdflatex manuscript.tex
	pdflatex manuscript.tex

# Convert and rename figures
figs:
	convert figures/Figure1.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure1.tiff
	convert figures/Figure2.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure2.tiff
	convert figures/Figure3.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure3.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_X_1_2_1_2-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure4.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure5.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-sgm_1-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure6.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-sgm_11-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure7.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-sgm_18-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure8.tiff
	convert figures/E-model-rifted-margin-grafics-estimated-model-alphas_X_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure9.tiff
	convert figures/E-model-rifted-margin-grafics-estimated-model-alphas_2_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure10.tiff
	convert figures/E-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure11.tiff
	convert figures/map-line-localization-test-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure12.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_X_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure13.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure14.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_7-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure15.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_17dot6-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure16.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_19-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure17.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_17-T07-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/Figure18.tiff

# Remove unused files
clean:
	rm *.aux *.bbl *.blg *.log *.pdf *.lof *.lot *.synctex.gz *.out

