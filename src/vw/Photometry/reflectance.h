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

#include "reconstruct.h"
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <math.h>

Vector3 computeNormalFrom3DPointsGeneral(Vector3 p1, Vector3 p2, Vector3 p3);
Vector3 computeNormalFrom3DPoints(Vector3 p1, Vector3 p2, Vector3 p3);
std ::vector<Vector3> ReadSunPosition(char *filename, int numEntries);
std ::vector<Vector3> ReadSpacecraftPosition(char *filename, int numEntries);
float computeReflectanceFromNormal(Vector3 sunPos, Vector3 xyz,  Vector3 normal);
float computeLambertianReflectanceFromNormal(Vector3 sunPos, Vector3 xyz,  Vector3 normal);
float computeLunarLambertianReflectanceFromNormal(Vector3 sunPos, Vector3 viewerPos, Vector3 xyz,  Vector3 normal, float B_0, float L);
float computeLunarLambertianReflectanceFromNormal(Vector3 sunPos, Vector3 viewPos, Vector3 xyz,  Vector3 normal);
float computeImageReflectance(std::string input_img_file, std::string DEM_file, std::string shadow_file, 
			      modelParams input_img_params, std::string output_img_file, GlobalParams globalParams);
float ComputeReflectance(Vector3 normal, Vector3 xyz, modelParams input_img_params, GlobalParams globalParams);
float computeImageReflectance(std::string input_img_file, std::string overlap_img_file, 
                              std::string DEM_file, 
                              std::string shadow_file, std::string overlap_shadow_file, 
                              modelParams input_img_params, modelParams overlap_img_params,
                              std::string output_img_file, GlobalParams globalParams);
