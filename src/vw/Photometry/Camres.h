// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Camres.h

#ifndef __VW_PHOTOMETRY_CAMRES_H__
#define __VW_PHOTOMETRY_CAMRES_H__

#include <string>
#include <vector>
#include <vw/Math/Vector.h>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {
#if 0
  int save_exposure_images(std::vector<std::string> output_files,
                           std::vector<std::string> input_files,
                           Vector<float> image_response, time_t mt_image_response);
  Vector<float> save_exposure_images(std::vector<std::string> output_files,
                                     std::vector<std::string> input_files,
                                     const char * input_file);
  int save_exposure_images(std::vector<std::string> output_files,
                           std::vector<std::string> input_files,
                           std::vector<std::string> camre_files);
#endif
}} // vw::photometry

#endif//__VW_PHOTOMETRY_CAMRES_H__
