INCPATH += -I../../src -I/opt/local/include
LIBS = ../../src/vw/FileIO/libvwFileIO.la
CXXFLAGS = -O3 -g -Wall

GENERATED_IMAGES = mural_clamp.jpg mural_crop.jpg			\
	mural_flip_horizontal.jpg mural_flip_vertical.jpg		\
	mural_rotate_180.jpg mural_rotate_90_ccw.jpg			\
	mural_rotate_90_cw.jpg mural_subsample.jpg mural_threshold.jpg	\
	mural_transpose.jpg

all: workbook.pdf

workbook.pdf: images
	pdflatex workbook.tex
	pdflatex workbook.tex

images: generate
	./generate

generate: generate.o
	../../libtool --mode=link $(CXX) -o $@ $^ $(LIBS)

.SUFFIXES: .cxx .c .o .cc

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

clean:
	rm -f *.o *~ \#* *.aux *.log $(GENERATED_IMAGES)

distclean: clean
	rm -rf generate workbook.pdf .libs

.PHONY: images
