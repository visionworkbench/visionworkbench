// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Outlier.h

#ifndef __VW_PHOTOMETRY_OUTLIER_H__
#define __VW_PHOTOMETRY_OUTLIER_H__

#include <vector>
#include <string>

namespace vw {
namespace photometry {

#if 0

  // Written by Taemin Kim
  float save_normal_images(std::vector<std::string> output_files,
                           std::vector<std::string> realexp_files,
                           std::vector<std::string> realrad_files,
                           char * exp_time_file);
  void weight_images(std::vector<std::string> output_files,
                     std::vector<std::string> input_files);

#endif

}} // no vw::photometry

#endif//__VW_PHOTOMETRY_OUTLIER_H__
