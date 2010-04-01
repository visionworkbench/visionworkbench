// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Shape.h

#ifndef __VW_PHOTOMETRY_SHAPE_H__
#define __VW_PHOTOMETRY_SHAPE_H__

#include <iostream>
#include <fstream>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>

#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  void InitDEM(ModelParams input_img_params,
               std::vector<ModelParams> overlap_img_params,
               GlobalParams globalParams);
  void ComputeSaveDEM(std::string curr_input_file, std::string prev_input_file,
                      std::string curr_output_file, std::string DEM_file,
                      ModelParams currModelParams,
                      ModelParams prevModelParams );

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_SHAPE_H__
