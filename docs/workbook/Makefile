INCPATH += -I../../src -I/opt/local/include -I/irg/packages/i386_darwin_gcc4.0/boost/include/
LIBS = ../../src/vw/libvw.la ../../src/vw/Mosaic/libvwMosaic.la -lboost_filesystem -L/irg/packages/i386_darwin_gcc4.0/boost/lib
CXXFLAGS = -O3 -g -Wall

GENERATED_IMAGES = images/mural_*.jpg images/hand.jpg images/lips.jpg	\
	images/hand-lips-blend.jpg mask.*.png images/Walker.qtree

all: workbook.pdf

workbook.pdf: images
	pdflatex workbook.tex
	pdflatex workbook.tex

images: generate
	rm -rf images/Walker.qtree
	./generate

generate: generate.o
	../../libtool --mode=link $(CXX) -o $@ $^ $(LIBS)

.SUFFIXES: .cxx .c .o .cc

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

clean:
	rm -rf *.o *~ \#* *.aux *.log $(GENERATED_IMAGES)

distclean: clean
	rm -rf generate workbook.toc workbook.pdf .libs

.PHONY: images
