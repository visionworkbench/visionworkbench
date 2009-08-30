TOOLS = ../../src/vw/tools

GENERATED_IMAGES = images/mural_*.jpg images/hand.jpg images/lips.jpg	\
	images/hand-lips-blend.jpg mask.*.png images/Walker.qtree       \
	images/pattern_*.jpg

all: workbook.pdf

workbook.pdf: images
	pdflatex workbook.tex
	pdflatex workbook.tex

images: $(TOOLS)/doc-generate
	rm -rf images/Walker.qtree
	$(TOOLS)/doc-generate

clean:
	rm -rf *.o *~ \#* *.aux *.log *.out $(GENERATED_IMAGES)

distclean: clean
	rm -rf generate workbook.toc workbook.pdf .libs

.PHONY: images clean distclean all
