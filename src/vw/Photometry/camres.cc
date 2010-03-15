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
#include "reconstruct.h"
#include "camres.h"
#include "misc.h"

int save_exposure_images(std::vector<std::string> output_files, std::vector<std::string> input_files, Vector<float> image_response, time_t mt_image_response) {
	int n=0;
	
	std::cout  << std::endl << " Constructing exposure images." << std::endl;
	for (unsigned i = 0; i < input_files.size(); ++i) {
		struct stat file_stat;
		if ( !stat(output_files[i].c_str(), &file_stat) ) 
			if ( difftime(file_stat.st_mtime, mt_image_response) > 0 ) {
				std::cout << " " << i << "th exposure image: " << output_files[i] << " is preserved." << std::endl;
				continue;
			}
		
		std::cout << " Constructing " << i << "th exposure image: " << output_files[i] << "." << std::endl;
		
		GeoReference geo;
		read_georeference(geo, input_files[i]);
		DiskImageView<PixelMask<PixelGray<uint8> > > image(input_files[i]);
		ImageView<PixelMask<PixelGray<float> > > tm_image(image.cols(), image.rows());
		for (unsigned x=0; x<image.cols(); ++x)
			for (unsigned y=0; y<image.rows(); ++y)
				if ( is_valid(image(x,y)) )
					tm_image(x,y) = image_response(image(x,y));
		
		tm_image = copy_mask(apply_mask(tm_image,0), image);
		write_georeferenced_image(output_files[i], tm_image, geo, TerminalProgressCallback("{Core}","Processing:"));
		n++;
	}
	return n;
}

Vector<float> save_exposure_images(std::vector<std::string> output_files, std::vector<std::string> input_files, const char * input_file) {
	float buffer[DYNAMIC_RANGE];
	static time_t mt_image_response;
	Vector<float> image_response(DYNAMIC_RANGE);
	
	struct stat file_stat;
	if ( stat(input_file, &file_stat) ) {
		std::cout  << std::endl << "Image response is linearly initialized." << std::endl;
		for (unsigned i=0; i<DYNAMIC_RANGE; ++i) {
			image_response(i)=i;
			buffer[i]=i;
		}
		save_binary_file(buffer,DYNAMIC_RANGE,input_file);	
		std::cout << "Image Response: " << image_response << std::endl;
		
		stat(input_file, &file_stat);
		save_exposure_images(output_files, input_files, image_response, file_stat.st_mtime);
	} else {
		load_binary_file(buffer,DYNAMIC_RANGE,input_file);	
		for (unsigned i=0; i<DYNAMIC_RANGE; ++i) image_response(i)=buffer[i];
		
		std::cout  << std::endl << "Reading Image response from the file " << input_file << "." << std::endl;
		if ( difftime(file_stat.st_mtime, mt_image_response) > 0 ) {
			mt_image_response = file_stat.st_mtime;
			std::cout << "Image response is modified at " << ctime(&file_stat.st_mtime);
			std::cout << "Image Response: " << image_response << std::endl;
			save_exposure_images(output_files, input_files, image_response, mt_image_response);
		} else {
			std::cout << "Image response is not modified since " << ctime(&mt_image_response);
		}
	}
	
	return image_response;
}

int save_exposure_images(std::vector<std::string> output_files, std::vector<std::string> input_files, std::vector<std::string> camre_files) {
	int n;
	GeoReference geo;
	float buffer[DYNAMIC_RANGE];
	Vector<float> image_response(DYNAMIC_RANGE);
	struct stat input_stat, output_stat, camre_stat;
	for (int i = 0; i < input_files.size(); ++i) {
		if ( !stat(camre_files[i].c_str(), &camre_stat) ) {
			if ( !stat(output_files[i].c_str(), &output_stat) ) {
				if ( difftime(output_stat.st_mtime, camre_stat.st_mtime) > 0 ) {
					std::cout << " " << i << "th exposure image: " << output_files[i] << " is preserved." << std::endl;			
					continue;		
				}
			}
			load_binary_file(buffer,DYNAMIC_RANGE,camre_files[i].c_str());	
			for (unsigned i=0; i<DYNAMIC_RANGE; ++i) image_response(i)=buffer[i];
		} else {
			std::cout  << std::endl << "Image response is linearly initialized." << std::endl;
			for (unsigned k=0; k<DYNAMIC_RANGE; ++k) {
				image_response(k)=k;
				buffer[k]=k;
			}
			save_binary_file(buffer,DYNAMIC_RANGE,camre_files[i].c_str());	
			std::cout << "Image Response: " << image_response << std::endl;
		}
		
		std::cout  << std::endl << "Reading Image response from the file " << camre_files[i] << "." << std::endl;
		std::cout << " " << i << "th exposure image: " << output_files[i] << " is generated." << std::endl;			
		read_georeference(geo, input_files[i]);
		DiskImageView<PixelMask<PixelGray<uint8> > > image(input_files[i]);
		ImageView<PixelMask<PixelGray<float> > > tm_image(image.cols(), image.rows());
		for (unsigned x=0; x<image.cols(); ++x)
			for (unsigned y=0; y<image.rows(); ++y)
				if ( is_valid(image(x,y)) )
					tm_image(x,y) = image_response(image(x,y));
		
		tm_image = copy_mask(apply_mask(tm_image,0), image);
		write_georeferenced_image(output_files[i], tm_image, geo, TerminalProgressCallback("{Core}","Processing:"));
		n++;
 	}
	return n;
}
