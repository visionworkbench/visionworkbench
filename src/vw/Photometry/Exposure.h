// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Exposure.h

#ifndef __VW_PHOTOMETRY_EXPOSURE_H__
#define __VW_PHOTOMETRY_EXPOSURE_H__

#include <string>
#include <vector>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  //computes the exposure time from image, albedo and DEM
  void ComputeExposure(ModelParams *currModelParams,
                       GlobalParams globalParams);

  //used for mosaicking, with no reflectance model
  void ComputeExposureAlbedo(ModelParams *currModelParams,
                             GlobalParams globalParams);
  /*
  void AppendExposureInfoToFile(std::string exposureFilename,
                                ModelParams currModelParams);
  std::vector<float> ReadExposureInfoFile(std::string exposureFilename,
                                          int numEntries);
  */
  void SaveExposureInfoToFile(ModelParams modelParams);
  void ReadExposureInfoFromFile(ModelParams *modelParams);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_EXPOSURE_H__
