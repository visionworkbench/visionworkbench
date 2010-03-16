#include <iostream>
#include <fstream>
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
// Written by Taemin Kim
float save_normal_images(std::vector<std::string> output_files, std::vector<std::string> realexp_files, 
						 std::vector<std::string> realrad_files, char * exp_time_file);
void weight_images(std::vector<std::string> output_files, std::vector<std::string> input_files);