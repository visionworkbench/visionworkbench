// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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

  int* ComputeImageCenterLine(std::string input_img_file,
                              int **r_maxDistArray);
  int* ComputeDEMCenterLine(std::string input_DEM_file, int noDataDEMVal,
                            int **r_maxDistArray);
  int* ComputeImageHorCenterLine(std::string input_img_file,
                                 int **r_maxDistArray);
  int* ComputeDEMHorCenterLine(std::string input_DEM_file, int noDataDEMVal,
                               int **r_maxDistArray);
  float ComputeLineWeights(Vector2 pix, int *centerLine,
                           int *maxDistArray);
  float ComputeLineWeightsV(Vector2 pix, int *centerLine,
                            int *maxDistArray);
  void SaveWeightsParamsToFile(struct ModelParams modelParams);
  void ReadWeightsParamsFromFile(struct ModelParams *modelParams);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_WEIGHTS_H__
