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

using namespace std;
using namespace vw;

#include <vw/Photometry/Reconstruct.h>

//subsamples a geo referenced tiff image by two 
void subsample_image(std::string output_file, std::string input_file);

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
                                   std::vector<Vector3> *xyzArray );


/// Erases a file suffix if one exists and returns the base string
std::vector<std::string> derivative_filenames(std::vector<std::string> input_files, std::string const& suffix);
void error_small_buffer(const char * file_name, unsigned entry_size);
template <class T>
void load_binary_file(T * ptr, size_t count, const char * filename);
template <class T>
void save_binary_file(T * ptr, size_t count, const char * filename);
Vector<uint8> load_inverse_weight(char * input_file, int num);
Vector2 get_minmax_values( std::string filename );
ImageView<PixelMask<PixelGray<float> > > interpolate_image(std::string input_file, std::string index_file, uint8 mask);
Vector<float> save_image_histograms(char * output_file, std::vector<std::string> output_files, std::vector<std::string> input_files,
                                                                        std::vector<std::string> index_files, std::vector<std::string> weight_files, char * weight_file);
std::vector<std::string> parse_command_arguments(int argc, char *argv[] );

Matrix<unsigned long> save_residual_histogram(char * output_file, std::vector<std::string> input_files,
                                                                                          std::vector<std::string> uintrad_files, Vector<float> exposure_times);
Matrix<unsigned long> save_residual_histogram(char * output_file, std::vector<std::string> input_files,
                                                                                          std::vector<std::string> uintrad_files, std::vector<std::string> weight_files, Vector<float> exposure_times);
Vector<float> save_logexp_histogram(char * output_file, std::vector<std::string> output_files, std::vector<std::string> input_files, std::vector<std::string> radiance_files,
                                                                        Vector<float> exposure_times, std::vector<std::string> index_files, std::vector<std::string> weight_files, char * weight_file);
float total_loglike(Vector<float> coeffs, Vector<float> image_histogram, Vector<float> image_response);
void convert_real2image(std::string input_file, std::string output_file, float ub, float lb = 0.0);
void convert_real2images(std::vector<std::string> input_files, std::vector<std::string> output_files, float ub, float lb = 0.0);
void convert_real2images(std::vector<std::string> input_files, std::vector<std::string> output_files, bool bNormalized = false);
