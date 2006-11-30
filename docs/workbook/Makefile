INCPATH += -I../../src -I/opt/local/include
LIBS = ../../src/vw/FileIO/libvwFileIO.la
CXXFLAGS = -O3 -g -Wall

GENERATED_IMAGES = images/mural_*.jpg images/hand.jpg images/lips.jpg	\
	images/hand-lips-blend.jpg mask.*.png

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
