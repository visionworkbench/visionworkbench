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