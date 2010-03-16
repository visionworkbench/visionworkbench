#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

//using namespace std;
//using namespace vw;
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions/gamma.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Index.h>

template <class ViewT>
std::vector<float> build_histogram(ImageViewBase<ViewT> const& view) {
	typedef typename ViewT::pixel_accessor pixel_accessor;
	std::vector<float> result(DYNAMIC_RANGE);
	
	int num_valid = 0;
	pixel_accessor row_acc = view.impl().origin();
	for( int32 row=0; row < view.impl().rows(); ++row ) {
		pixel_accessor col_acc = row_acc;
		for( int32 col=0; col < view.impl().cols(); ++col ) {
			if ( is_valid(*col_acc) ) {
				result[ (*col_acc)[0] ] += 1.0;
				++num_valid;
			}
			col_acc.next_col();
		}
		row_acc.next_row();
	}
	
	// Renormalize the histogram values
	for (int32 i=0; i< result.size(); ++i) 
		result[i] /= num_valid;
	
	return result;
}

void display_histogram(std::string input_file) {
	DiskImageView<PixelMask<PixelGray<uint8> > > image(input_file);
	std::vector<float> hist = build_histogram(image);
	std::cout << "Histogram of " << input_file << " ";
	for (unsigned i=0; i < hist.size(); ++i) {
		if (hist[i]) std::cout << i << ":" << hist[i] << " ";
	}
	std::cout << std::endl;
}

void index_image(std::string index_file, std::string input_file, uint8 mask) {
	GeoReference geo1, geo2;
	read_georeference(geo1, index_file);
	read_georeference(geo2, input_file);
	DiskImageView<PixelMask<PixelGray<uint8> > > image(input_file);
	DiskImageView<PixelMask<PixelGray<uint8> > > index(index_file);
	ImageView<PixelMask<PixelGray<uint8> > > tm_index = index;
	
	for (unsigned x=0; x<index.cols(); ++x)
		for (unsigned y=0; y<index.rows(); ++y)
			if ( is_valid(index(x,y)) ) {
				Vector2 subpix = geo2.lonlat_to_pixel(geo1.pixel_to_lonlat(Vector2(x,y)));
				int i=int(subpix[0]);
				int j=int(subpix[1]);
				if ( i < 0 || j < 0 || i > image.cols()-2 || j > image.rows()-2 ) continue;
				if ( is_valid(image(i,j)) && is_valid(image(i+1,j)) && is_valid(image(i,j+1)) && is_valid(image(i+1,j+1)) ) tm_index(x,y) += mask;
			}
	write_georeferenced_image(index_file, tm_index, geo1, TerminalProgressCallback("{Core}","Processing:"));
}

void index_images(std::vector<std::string> index_files, std::vector<std::string> input_files) {
	for (unsigned i = 0; i < input_files.size(); ++i) {
		struct stat index_stat;
		if ( stat(index_files[i].c_str(), &index_stat) ) { 
			GeoReference geo1;
			read_georeference(geo1, input_files[i]);
			DiskImageView<PixelMask<PixelGray<uint8> > > image(input_files[i]);
			ImageView<PixelMask<PixelGray<uint8> > > index(image.cols(), image.rows());
			index = copy_mask(apply_mask(index,0), image);
			
			std::cout << " " << i << "th index image: " << index_files[i] << " is created." << std::endl;
			write_georeferenced_image(index_files[i], index, geo1, TerminalProgressCallback("{Core}","Processing:"));
			
			int jmin = (i>4) ? i-4: 0;
			for (int j = jmin; j<i; ++j) {
				std::cout << "(" << i << "," << j << ") " << index_files[i].c_str() << " = " << input_files[i].c_str() << " & " << input_files[j].c_str() << std::endl;
				index_image(index_files[i], input_files[j], 0x01 << j-i+4);
			}
			
			int jmax = ((i+5)<input_files.size()) ? i+5: input_files.size();
			for (int j = i+1; j<jmax; ++j) {
				std::cout << "(" << i << "," << j << ") " << index_files[i].c_str() << " = " << input_files[i].c_str() << " & " << input_files[j].c_str() << std::endl;
				index_image(index_files[i], input_files[j], 0x10 << j-i-1);
			}
			display_histogram(index_files[i]);
		}
	}
}
