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
	convert figures/rifted-margin-model-color.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/rifted-margin-model-color.tiff
	convert figures/interpretation-model-color.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/interpretation-model-color.tiff
	convert figures/flow_diagram.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/flow_diagram.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_X_1_2_1_2-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/B_X_1_2_1_2.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/B_3_1_2_1_2.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-sgm_1-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/B_3_1_2_1_2_sgm1.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-sgm_11-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/B_3_1_2_1_2_sgm11.tiff
	convert figures/B-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-sgm_18-2pts-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/B_3_1_2_1_2_sgm18.tiff
	convert figures/E-model-rifted-margin-grafics-estimated-model-alphas_X_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/E_X_1_2_1_2.tiff
	convert figures/E-model-rifted-margin-grafics-estimated-model-alphas_2_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/E_2_1_2_1_2.tiff
	convert figures/E-model-rifted-margin-grafics-estimated-model-alphas_3_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/E_3_1_2_1_2.tiff
	convert figures/map-line-localization-test-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/line_localization.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_X_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/pelotas_X_1_2_1_2.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/pelotas_2_1_2_1_2.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_7-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/pelotas_2_1_2_1_2_sgm7.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_17dot6-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/pelotas_2_1_2_1_2_sgm17dot6.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_19-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/pelotas_2_1_2_1_2_sgm19.tiff
	convert figures/pelotas-profile-grafics-estimated-model-alphas_2_1_2_1_2-sgm_17-T07-dpi300.png -colorspace CMYK -units PixelsPerInch -resample 400 -resize 1732 figures/pelotas_2_1_2_1_2_sgm17_alternative.tiff

# Remove unused files
clean:
	rm *.aux *.bbl *.blg *.log *.pdf *.lof *.lot *.synctex.gz *.out

