// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ReconstructError.h

#ifndef __VW_PHOTOMETRY_RECONSTRUCTERROR_H__
#define __VW_PHOTOMETRY_RECONSTRUCTERROR_H__

#include <string>
#include <vector>
//#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  float ComputeError(float intensity, float T,
                          float albedo, float reflectance);
                          //Vector3 /*xyz*/, Vector3 /*xyz_prior*/)

  //reconstruction error functions
  void ComputeReconstructionErrorMap(ModelParams input_img_params,
                                     std::vector<ModelParams> overlap_img_params,
                                     GlobalParams globalParams,
                                     float *avgError, int *totalNumSamples);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_RECONSTRUCTERROR_H__
