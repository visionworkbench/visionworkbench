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
//computes the exposure time from image, albedo and DEM
void ComputeExposure(std::string curr_input_file,
                     std::string curr_albedo_file,
                     std::string DEM_file,
                     modelParams *currModelParams,
                     GlobalParams globalParams);

/*
//used for mosaicking, with no reflectance model
void ComputeExposure(std::string curr_input_file,
                     std::string curr_albedo_file,
                     modelParams *currModelParams,
                     GlobalParams globalParams);
*/

//used for mosaicking, with no reflectance model
void ComputeExposure(modelParams *currModelParams,
                     GlobalParams globalParams);

void AppendExposureInfoToFile(string exposureFilename, string currInputFile, modelParams currModelParams);
vector<float> ReadExposureInfoFile(string exposureFilename, int numEntries);

//old function, no longer used
//void ComputeExposureTime(std::vector<Vector4> const& samples, modelParams prevModelParams, modelParams *currModelParams,
//                         std::vector<Vector3> normalArray, std::vector<Vector3> xyzArray, float *error);

// Written by Taemin Kim
Vector<float> save_exposure_times(char * output_file, std::vector<std::string> radiance_files, std::vector<std::string> response_files,
                                                                  std::vector<std::string> index_files, char * weight_file);
// Written by Taemin Kim
// Preserve output_files if output_files is modified later than input_files and exposure_times.
Vector<float> save_exposure_times(char * output_file, std::vector<std::string> radiance_files, std::vector<std::string> response_files,
                                                                  std::vector<std::string> index_files, std::vector<std::string> weight_files, char * weight_file);

Vector<float> load_exposure_times(char * input_file, int num);



