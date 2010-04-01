// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Camres.h

#ifndef __VW_PHOTOMETRY_CAMRES_H__
#define __VW_PHOTOMETRY_CAMRES_H__

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

int save_exposure_images(std::vector<std::string> output_files, std::vector<std::string> input_files,
                                                 Vector<float> image_response, time_t mt_image_response);
Vector<float> save_exposure_images(std::vector<std::string> output_files, std::vector<std::string> input_files,
                                                                   const char * input_file);
int save_exposure_images(std::vector<std::string> output_files, std::vector<std::string> input_files,
                                                 std::vector<std::string> camre_files);

#endif//__VW_PHOTOMETRY_CAMRES_H__
