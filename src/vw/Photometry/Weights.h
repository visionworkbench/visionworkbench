// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Weights.h

#ifndef __VW_PHOTOMETRY_WEIGHTS_H__
#define __VW_PHOTOMETRY_WEIGHTS_H__

#include <vw/Math/Vector.h>
#include <string>

namespace vw {
namespace photometry {

  int* ComputeImageHCenterLine(std::string input_img_file,
                              int **r_hMaxDistArray);
  int* ComputeDEMHCenterLine(std::string input_DEM_file, int noDataDEMVal,
                            int **r_hMaxDistArray);
  int* ComputeImageVCenterLine(std::string input_img_file,
                                 int **r_hMaxDistArray);
  int* ComputeDEMVCenterLine(std::string input_DEM_file, int noDataDEMVal,
                               int **r_hMaxDistArray);
  float ComputeLineWeightsH(Vector2 pix, int *hCenterLine,
                           int *hMaxDistArray);
  float ComputeLineWeightsV(Vector2 pix, int *hCenterLine,
                            int *hMaxDistArray);
  float ComputeLineWeightsHV(Vector2 pix,
                             int *hCenterLine, int *hMaxDistArray,
                             int *vCenterLine, int *vMaxDistArray
                             );
  void SaveWeightsParamsToFile(struct ModelParams modelParams);
  void ReadWeightsParamsFromFile(struct ModelParams *modelParams);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_WEIGHTS_H__
