// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

//#include "shape_from_shading.h"
//#include "shape.h"
//#include "albedo.h"
//#include "exposure.h"
//#include "reflectance.h"
//#include "reconstruct.h"
//#include "shadow.h"
//#include "index.h"
//#include <math.h>
//#include <cv.h>
//#include <highgui.h>

#include "reconstruct.h"

//Vector2 ComputeImageCenter2D(std::string input_img_file, float *maxDist);
//float ComputeWeights(Vector2 pix, Vector2 C, float maxDistance);
int* ComputeImageCenterLine(std::string input_img_file, int **r_maxDistArray);
int* ComputeDEMCenterLine(std::string input_DEM_file, int **r_maxDistArray);
int* ComputeImageHorCenterLine(std::string input_img_file, int **r_maxDistArray);
int* ComputeDEMHorCenterLine(std::string input_DEM_file, int **r_maxDistArray);
float ComputeLineWeights(Vector2 pix, int *centerLine, int *maxDistArray);
float ComputeLineWeightsV(Vector2 pix, int *centerLine, int *maxDistArray);


