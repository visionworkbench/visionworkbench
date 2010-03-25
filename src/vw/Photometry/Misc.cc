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

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Misc.h>


#include <string>
#include <fstream>
#include <complex>
#include <vector>
#include <limits>
#include <time.h>
#include <unistd.h>
#include <proj_api.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po=boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

//subsamples a geo referenced tiff image by two 
void subsample_image(std::string output_file, std::string input_file) {
	GeoReference geo;
	read_georeference(geo, input_file);
        DiskImageView<PixelMask<PixelGray<uint8> > >  image(input_file); 
	int cols = (image.cols()+1)/2, rows = (image.rows()+1)/2;
	ImageView<PixelMask<PixelGray<uint8> > > tm_image(cols, rows);

      
        ImageViewRef<PixelMask<PixelGray<uint8> > >  interp = interpolate(edge_extend(image.impl(), 
                                                                         ConstantEdgeExtension()), 
                                                                         BilinearInterpolation());
       
	int x, y;

	for (x=0; x<cols; ++x){
	  for (y=0; y<rows; ++y){
	      //if ( is_valid(image(2*x,2*y)) || is_valid(image(2*x+1,2*y)) || is_valid(image(2*x,2*y+1)) || is_valid(image(2*x+1,2*y+1)) ) {
	      //if (is_valid(image(2*x,2*y)) && is_valid(image(2*x+1,2*y)) && is_valid(image(2*x,2*y+1)) && is_valid(image(2*x+1,2*y+1)) ) {
	    if ( is_valid(interp(2*x+0.5, 2*y+0.5)) ){
		  tm_image(x,y) = interp(2*x+0.5, 2*y+0.5);
	      } 
              else{
		  tm_image(x,y).invalidate();
	      }
            }
	}
	
	Matrix<double> H = geo.transform();
	H(0,0) *= 2;
	H(1,1) *= 2;
	geo.set_transform(H);
	
	write_georeferenced_image(output_file, tm_image, geo, TerminalProgressCallback("{Core}","Processing:"));
}


// Given two images and two georeferences, this function picks a set
// of matching pixel samples between the two images.  It rejects
// pixels that are not valid, and it should probably also reject
// pixels that are near saturation (though it does not yet!).
template <class ViewT>
std::vector<Vector4> sample_images(ImageViewBase<ViewT> const& image1, 
                                   ImageViewBase<ViewT> const& image2,
                                   GeoReference const& geo1,
                                   GeoReference const& geo2,
                                   int num_samples,
                                   std::string const& DEM_file, 
                                   std::vector<Vector3> *normalArray,
                                   std::vector<Vector3> *xyzArray ) {
  int sample = 0;
  int numtries = 0;
  std::vector<Vector4> result;
  
  // Random numbers
  srandom((unsigned int) clock());
  float rand_max = pow(2.0,31)-1;
  
  ImageViewRef<typename ViewT::pixel_type> interp_image1 = interpolate(edge_extend(image1.impl(), 
                                                                       ConstantEdgeExtension()), 
                                                                       BilinearInterpolation());
  ImageViewRef<typename ViewT::pixel_type> interp_image2 = interpolate(edge_extend(image2.impl(), 
                                                                       ConstantEdgeExtension()), 
                                                                       BilinearInterpolation());

  // This block of code samples the images comprehensively, adding a
  // sample pair for every valid pixel in interp_image1.
  // for (unsigned j=0; j < interp_image1.rows(); ++j) {
  //   for (unsigned i=0; i < interp_image1.cols(); ++i) {
  //     Vector2 sample_pix1(i,j);
  //     Vector2 sample_pix2 = geo2.lonlat_to_pixel(geo1.pixel_to_lonlat(sample_pix1));
      
  //     // Check to see whether these pixels are valid
  //     typename ViewT::pixel_type pix1 = interp_image1(sample_pix1[0], sample_pix1[1]);
  //     typename ViewT::pixel_type pix2 = interp_image2(sample_pix2[0], sample_pix2[1]);
  //     if ( is_valid(pix1) && is_valid(pix2) &&
  //          pix1[0] > 10 && pix1[0] < 245 &&
  //          pix2[0] > 10 && pix2[0] < 245 ) {
  //       result.push_back(Vector2(pix1[0],pix2[0]));
  //       ++sample;
  //       //        std::cout << result[result.size()-1][0] << " " << result[result.size()-1][1] << "\n";
  //     }
  //     ++numtries;
  //   }
  // }
  

  //added by Ara to support DEMs - START
  DiskImageView<PixelGray<float> >  dem_image(DEM_file); 
  GeoReference GR;
  read_georeference(GR, DEM_file);
  //added by Ara to support DEMs - END

  // This block of code samples the images randomly, gathering up to
  // num_samples samples from the images. 
  while (sample < num_samples && numtries < num_samples*10) {
    
    Vector2 sample_pix1(float(random())/rand_max * image1.impl().cols(),
                        float(random())/rand_max * image1.impl().rows());
    Vector2 sample_pix2 = geo2.lonlat_to_pixel(geo1.pixel_to_lonlat(sample_pix1));

    Vector2 sample_pix_dem = GR.lonlat_to_pixel(geo1.pixel_to_lonlat(sample_pix1));

    Vector2 lonlat = geo1.pixel_to_lonlat(sample_pix1);
   
    // Check to see whether these pixels are valid
    typename ViewT::pixel_type pix1 = interp_image1(sample_pix1[0], sample_pix1[1]);
    typename ViewT::pixel_type pix2 = interp_image2(sample_pix2[0], sample_pix2[1]);
    if ( is_valid(pix1) && is_valid(pix2) ) {

       //result.push_back(Vector4(pix1[0],pix2[0],lonlat[0],lonlat[1]));
       
       int x = (int)sample_pix_dem[0];
       int y = (int)sample_pix_dem[1];

       if (x < 0){
           x = 0;
       }
       if (x > dem_image.cols()-1){
	   x = dem_image.cols()-1; 
       }
       if (y < 0){
           y = 0;
       }
       if (y > dem_image.rows()-1){
	   y = dem_image.rows()-1; 
       }
      
       Vector3 longlat3(lonlat(0),lonlat(1),(dem_image)(x, y));
       Vector3 xyz = geo1.datum().geodetic_to_cartesian(longlat3);
       
       Vector2 sample_pix_dem_left;
       sample_pix_dem_left(0) = x-1;
       if (sample_pix_dem_left(0) < 0){
	  sample_pix_dem_left(0) = 0;
	  //break;
       }
       sample_pix_dem_left(1) = y;
       lonlat = GR.pixel_to_lonlat(sample_pix_dem_left);
    
       Vector3 longlat3_left(lonlat(0),lonlat(1),(dem_image)(sample_pix_dem_left(0), sample_pix_dem_left(1)));
       Vector3 xyz_left = geo1.datum().geodetic_to_cartesian(longlat3_left);
       
       Vector2 sample_pix_dem_top;
       sample_pix_dem_top(0) = x;
       sample_pix_dem_top(1) = y-1;
       if (sample_pix_dem_top(1) < 0){
	 sample_pix_dem_top(1) = 0;
	 //break;
       }

       lonlat = GR.pixel_to_lonlat(sample_pix_dem_top);     
       Vector3 longlat3_top(lonlat(0),lonlat(1),(dem_image)(sample_pix_dem_top(0), sample_pix_dem_top(1)));
       Vector3 xyz_top = geo1.datum().geodetic_to_cartesian(longlat3_top);
       
       

       Vector3 normal = computeNormalFrom3DPoints(xyz, xyz_left, xyz_top);

       //printf("normal:%f %f %f\n", normal(0), normal(1), normal(2));
       //printf("%f %f %f\n", loc_longlat3(0), loc_longlat3(1), loc_longlat3(2));
       
       result.push_back(Vector4(pix1[0],pix2[0],lonlat[0],lonlat[1]));
        
       normalArray->push_back(normal);
       xyzArray->push_back(xyz);
       
       ++sample;
       //      std::cout << result[result.size()-1][0] << " " << result[result.size()-1][1] << "\n";
    } 
    ++numtries;
  }
  return result;
}

/// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
        std::string result = filename;
        int index = result.rfind(".");
        if (index != -1)
                result.erase(index, result.size());
        return result;
}

/// Erases a file suffix if one exists and returns the base string
std::vector<std::string> derivative_filenames(std::vector<std::string> input_files, std::string const& suffix) {
        std::vector<std::string> output_files(input_files.size());
        for (unsigned i=0; i<input_files.size(); ++i)
                output_files[i]=prefix_from_filename(input_files[i]) + suffix;
        return output_files;
}

void error_small_buffer(const char * file_name, unsigned entry_size) {
        std::cout << "The size of buffer is limited by " << SIZE_OF_BUFFER << "." << std::endl;
        std::cout << "The number of entry is " << entry_size << ", which is bigger than buffer size." << std::endl;
        std::cout << "The entries from / in the file, " << file_name << " , may contain garbages!!!" << std::endl;
}

template <class T>
void load_binary_file(T * ptr, size_t count, const char * filename) {
        FILE * fp = fopen(filename, "rb");
        if (fp != NULL) {
                fread(ptr, sizeof(T), count, fp);
                fclose(fp);
                if (count > SIZE_OF_BUFFER) error_small_buffer(filename, count);
        } else {
                std::cout << "Error to read " << filename << std::endl;
                exit(1);
        }
}

template <class T>
void save_binary_file(T * ptr, size_t count, const char * filename) {
        FILE * fp = fopen(filename, "wb");
        if (fp != NULL) {
                fwrite(ptr, sizeof(T), count, fp);
                fclose(fp);
                if (count > SIZE_OF_BUFFER) error_small_buffer(filename, count);
        } else {
                std::cout << "Error to write " << filename << std::endl;
                exit(1);
        }
}

Vector<uint8> load_inverse_weight(char * input_file, int num) {
        Vector<uint8> inverse_weight(num);
        float buffer[SIZE_OF_BUFFER];
        //FILE *fp;

        struct stat file_stat;
        if ( stat(input_file, &file_stat) ) {
                std::cout << "Generate inverse weight to  " << input_file << std::endl;
                for (int i = 0; i < num; ++i) {
                        inverse_weight(i) = 1;
                        for (unsigned j = 0; j < 8; ++j)
                                inverse_weight(i) += (0x01 & (i>>j));
                        buffer[i] = inverse_weight(i);
                        //                        buffer[i] = 1;
                }
                Vector<int> tmp = inverse_weight;
                std::cout << tmp << std::endl;
                save_binary_file(buffer,num,input_file);
        } else {
                std::cout << "Reading inverse weight from " << input_file << std::endl;
                load_binary_file(buffer, num, input_file);
                for (int i = 0; i < num; ++i) inverse_weight(i) = buffer[i];
                Vector<int> tmp = inverse_weight;
                std::cout << tmp << std::endl;
        }

        return inverse_weight;
}

Vector2 get_minmax_values( std::string filename ) {
        std::string cache_filename = prefix_from_filename(filename) + ".minmax";
        uint8 lo, hi;

        if ( fs::exists( cache_filename ) ) {
                std::ifstream input_file(cache_filename.c_str());
                int a, b;
                input_file >> a >> b;
                input_file.close();
                lo = a;
                hi = b;
        } else {
                DiskImageView<PixelMask<PixelGray<uint8> > > image(filename);
                min_max_channel_values(image, lo, hi);
                std::ofstream output_file(cache_filename.c_str());
                output_file << int(lo) << " " << int(hi);
                output_file.close();
        }
        return Vector2(lo, hi);
}

ImageView<PixelMask<PixelGray<float> > > interpolate_image(std::string input_file, std::string index_file, uint8 mask) {
        GeoReference geo1, geo2;
        read_georeference(geo1, index_file);
        read_georeference(geo2, input_file);
        DiskImageView<PixelMask<PixelGray<uint8> > > index(index_file);
        DiskImageView<PixelMask<PixelGray<float> > > image(input_file);
        ImageView<PixelMask<PixelGray<float> > > tm_image(index.cols(), index.rows());
        ImageViewRef<PixelMask<PixelGray<float> > > interp = interpolate(image, BilinearInterpolation());

        for (int x=0; x<(int)index.cols(); ++x)
	  for (int y=0; y<(int)index.rows(); ++y)
                        if ( is_valid(index(x,y)) && ( uint8(index(x,y)) & mask )) {
                                Vector2 subpix = geo2.lonlat_to_pixel(geo1.pixel_to_lonlat(Vector2(x,y)));
                                tm_image(x,y) = interp(subpix[0], subpix[1]);
                        }

        tm_image = copy_mask(apply_mask(tm_image,0), index);
        return tm_image;
}

Vector<float> save_image_histograms(char * output_file, std::vector<std::string> output_files, std::vector<std::string> input_files,
                                                                        std::vector<std::string> index_files, std::vector<std::string> weight_files, char * weight_file) {
        float buffer[DYNAMIC_RANGE];
        Vector<float> images_histogram(DYNAMIC_RANGE);
        Vector<uint8> inverse_weight = load_inverse_weight(weight_file, 256);

        for (unsigned i = 0; i < input_files.size(); ++i) {
                Vector<float> image_histogram(DYNAMIC_RANGE);
                DiskImageView<PixelMask<PixelGray<uint8> > > image(input_files[i]);
                DiskImageView<PixelMask<PixelGray<uint8> > > index(index_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > weight(weight_files[i]);
                for (int x=0; x<(int)image.cols(); ++x)
		  for (int y=0; y<(int)image.rows(); ++y)
                                if ( is_valid(image(x,y)) ) image_histogram(image(x,y)) += weight(x,y)/inverse_weight(index(x,y));

                for (unsigned k = 0; k < DYNAMIC_RANGE; ++k) buffer[k] = image_histogram(k);
                save_binary_file(buffer,DYNAMIC_RANGE,output_files[i].c_str());
                std::cout << "Histogram of " << input_files[i] << ": " << image_histogram << std::endl;

                images_histogram += image_histogram;
        }

        for (unsigned i = 0; i < DYNAMIC_RANGE; ++i) buffer[i] = images_histogram(i);
        save_binary_file(buffer,DYNAMIC_RANGE,output_file);
        std::cout << "Histogram of " << output_file << ": " << images_histogram << std::endl;

        return images_histogram;
}

// Create the output, index, and radiance file names

std::vector<std::string> parse_command_arguments(int argc, char *argv[] ) {
        int num_matches;
        std::vector<std::string> input_files;

        po::options_description general_options("Options");
        general_options.add_options()
        ("help", "Display this help message")
        ("num-matches,m", po::value<int>(&num_matches)->default_value(1000), "Number of points to match for linear regression.");

        po::options_description hidden_options("");
        hidden_options.add_options()
        ("input-files", po::value<std::vector<std::string> >(&input_files));

        po::options_description options("Allowed Options");
        options.add(general_options).add(hidden_options);

        po::positional_options_description p;
        p.add("input-files", -1);

        po::variables_map vm;
        po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
        po::notify( vm );

        std::ostringstream usage;
        usage << "Description: tonematches several images" << std::endl << std::endl;
        usage << "Usage: histeq [options] <filename1> <filename2> ..." << std::endl << std::endl;
        usage << general_options << std::endl;

        if( vm.count("help") ) {
                std::cerr << usage.str() << std::endl;
                exit(1);
        }

        if( vm.count("input-files")<1 ) {
                std::cerr << "Error: Must specify at least one input file!" << std::endl << std::endl;
                std::cerr << usage.str();
                exit(1);
        }

        return input_files;
}

Matrix<unsigned long> save_residual_histogram(char * output_file, std::vector<std::string> input_files,
                                                                                          std::vector<std::string> uintrad_files, Vector<float> exposure_times) {
        Matrix<unsigned long> histo(DYNAMIC_RANGE,DYNAMIC_RANGE);
        uint8 radiance_estimate;

        std::cout << "Compute Residual of Images: ";
        for (unsigned i = 0; i < input_files.size(); ++i) {
                DiskImageView<PixelMask<PixelGray<float> > > image(input_files[i]);
                DiskImageView<PixelMask<PixelGray<uint8> > > radiance(uintrad_files[i]);
                for (int x=0; x<(int)image.cols(); ++x)
		  for (int y=0; y<(int)image.rows(); ++y)
                                if ( is_valid(image(x,y)) ) {
                                        radiance_estimate = uint8(image(x,y)/exposure_times[i]);
                                        ++histo(radiance(x,y),radiance_estimate);
                                }
                std::cout << i << " ";
        }
        std::cout << std::endl;

        float buffer[DYNAMIC_RANGE][DYNAMIC_RANGE];
        for (unsigned i = 0; i < DYNAMIC_RANGE; ++i)
                for (unsigned j = 0; j < DYNAMIC_RANGE; ++j)
                        buffer[i][j] = histo(i,j);
        save_binary_file(buffer,DYNAMIC_RANGE*DYNAMIC_RANGE,output_file);

        return histo;
}

Matrix<unsigned long> save_residual_histogram(char * output_file, std::vector<std::string> input_files,
                                                                                          std::vector<std::string> uintrad_files, std::vector<std::string> weight_files, Vector<float> exposure_times) {
        Matrix<unsigned long> histo(DYNAMIC_RANGE,DYNAMIC_RANGE);
        uint8 radiance_estimate;

        std::cout << "Compute Residual of Images: ";
        for (unsigned i = 0; i < input_files.size(); ++i) {
                DiskImageView<PixelMask<PixelGray<float> > > image(input_files[i]);
                DiskImageView<PixelMask<PixelGray<uint8> > > radiance(uintrad_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > weight(weight_files[i]);
                for (int x=0; x<(int)image.cols(); ++x)
		     for (int y=0; y<(int)image.rows(); ++y)
                                if ( is_valid(image(x,y)) ) {
                                        radiance_estimate = uint8(image(x,y)/exposure_times[i]);
                                        histo(radiance(x,y),radiance_estimate) += weight(x,y);
                                }
                std::cout << i << " ";
        }
        std::cout << std::endl;

        float buffer[DYNAMIC_RANGE][DYNAMIC_RANGE];
        for (unsigned i = 0; i < DYNAMIC_RANGE; ++i)
                for (unsigned j = 0; j < DYNAMIC_RANGE; ++j)
                        buffer[i][j] = histo(i,j);
        save_binary_file(buffer,DYNAMIC_RANGE*DYNAMIC_RANGE,output_file);

        return histo;
}

Vector<float> save_logexp_histogram(char * output_file, std::vector<std::string> output_files, std::vector<std::string> input_files, std::vector<std::string> radiance_files,
                                                                        Vector<float> exposure_times, std::vector<std::string> index_files, std::vector<std::string> weight_files, char * weight_file) {
        float sensor_exposure;
        float buffer[DYNAMIC_RANGE];
        Vector<float> images_histogram(DYNAMIC_RANGE+1);
        Vector<uint8> inverse_weight = load_inverse_weight(weight_file, 256);

        for (unsigned i = 0; i < input_files.size(); ++i) {
                Vector<float> image_histogram(DYNAMIC_RANGE+1);
                DiskImageView<PixelMask<PixelGray<uint8> > > image(input_files[i]);
                DiskImageView<PixelMask<PixelGray<uint8> > > index(index_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > radiance(radiance_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > weight(weight_files[i]);
                for (int x=0; x<(int)image.cols(); ++x)
		  for (int y=0; y<(int)image.rows(); ++y)
                                if ( is_valid(image(x,y)) ) {
                                        sensor_exposure = radiance(x,y)*exposure_times[i];
                                        image_histogram(image(x,y)) += weight(x,y)*log(sensor_exposure)/inverse_weight(index(x,y));
                                        image_histogram(DYNAMIC_RANGE) += weight(x,y)*sensor_exposure/inverse_weight(index(x,y));
                                }
                for (unsigned k = 0; k < DYNAMIC_RANGE; ++k) buffer[k] = image_histogram(k);
                save_binary_file(buffer,DYNAMIC_RANGE,output_files[i].c_str());
                std::cout << "Logarithmic Exposure Histogram of " << input_files[i] << ": " << image_histogram << std::endl;

                images_histogram += image_histogram;
        }

        for (unsigned i = 0; i < DYNAMIC_RANGE; ++i) buffer[i] = images_histogram(i);
        save_binary_file(buffer,DYNAMIC_RANGE,output_file);
        std::cout << "Logarithmic Exposure Histogram of " << output_file << ": " << images_histogram << std::endl;

        return images_histogram;
}

float total_loglike(Vector<float> coeffs, Vector<float> image_histogram, Vector<float> image_response) {
        float sum_loglike = -coeffs(DYNAMIC_RANGE);

        std::cout << "init  " << sum_loglike << std::endl;
        for (unsigned i = 0; i < DYNAMIC_RANGE; ++i) {
                sum_loglike += coeffs(i)*image_response(i) - image_histogram(i)*lgamma(image_response(i)+1);
                std::cout << i << " " << sum_loglike << " " << coeffs(i) << " " << image_response(i) << " " <<  image_histogram(i) << " " << lgamma(image_response(i)+1) << std::endl;
        }
        return sum_loglike;
}

void convert_real2image(std::string input_file, std::string output_file, float ub, float lb) {
        GeoReference geo;
        read_georeference(geo, input_file);
        DiskImageView<PixelMask<PixelGray<float> > > image(input_file);

        if (ub == DYNAMIC_RANGE-1 && lb == 0.0) {
                write_georeferenced_image(output_file, channel_cast<uint8>(clamp(image,0.0,255.0)), geo, TerminalProgressCallback("{Core}","Processing:"));
        } else {
                ImageView<PixelMask<PixelGray<uint8> > > tm_image = (image-lb)*(DYNAMIC_RANGE-1)/(ub-lb);

                std::cout << input_file << " in [ " << lb << ", " << ub << "] => " << output_file << "." << std::endl;
                write_georeferenced_image(output_file, channel_cast<uint8>(clamp(tm_image, 0.0,255.0)), geo, TerminalProgressCallback("{Core}","Processing:"));
        }
}

void convert_real2images(std::vector<std::string> input_files, std::vector<std::string> output_files, float ub, float lb) {
        for (unsigned i = 0; i < input_files.size(); ++i) {
                std::cout << " " << i << "th image: " ;
                convert_real2image(input_files[i], output_files[i], ub, lb);
        }
}

void convert_real2images(std::vector<std::string> input_files, std::vector<std::string> output_files, bool bNormalized) {
        float lb = DYNAMIC_RANGE-1, ub = 0, lo, hi;
        if (bNormalized) {
                for (unsigned i = 0; i < input_files.size(); ++i) {
                        DiskImageView<PixelMask<PixelGray<float> > > image(input_files[i]);
                        min_max_channel_values(image, lo, hi);
                        if (lb > lo) lb = lo;
                        if (ub < hi) ub = hi;
                        std::cout << "Range of the " << i << "th image " << input_files[i] << " : (" << lo << "," << hi << ")" << std::endl;
                }

                std::cout << "Range of all image: (" << lb << "," << ub << ")" << std::endl;
        } else {
                lb =0;
                ub = DYNAMIC_RANGE-1;
        }

        for (unsigned i = 0; i < input_files.size(); ++i) {
                std::cout << " " << i << "th image: " ;
                convert_real2image(input_files[i], output_files[i], ub, lb);
        }
}



