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
/*
void ComputeExposure(std::string curr_input_file,
                     std::string curr_albedo_file,
                     std::string DEM_file,
                     ModelParams *currModelParams,
                     GlobalParams globalParams);
*/
void ComputeExposure(ModelParams *currModelParams,
                     GlobalParams globalParams);

//used for mosaicking, with no reflectance model
void ComputeExposureAlbedo(ModelParams *currModelParams,
                           GlobalParams globalParams);


void AppendExposureInfoToFile(string exposureFilename, ModelParams currModelParams);
vector<float> ReadExposureInfoFile(string exposureFilename, int numEntries);



