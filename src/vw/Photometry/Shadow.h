// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Shadow.h

#ifndef __VW_PHOTOMETRY_SHADOW_H__
#define __VW_PHOTOMETRY_SHADOW_H__

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  void ComputeSaveShadowMap( ModelParams input_img_params,
                             GlobalParams globalParams);
  void AddShadows(std::string input_img_file,
                  std::string output_img_file,
                  std::string shadow_file);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_SHADOW_H__
