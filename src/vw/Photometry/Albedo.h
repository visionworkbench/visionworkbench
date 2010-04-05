// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Albedo.h

#ifndef __VW_PHOTOMETRY_ALBEDO_H__
#define __VW_PHOTOMETRY_ALBEDO_H__

#include <string>
#include <vector>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  //image mosaic functions
  void InitImageMosaic(ModelParams input_img_params,
                       std::vector<ModelParams> overlap_img_params,
                       GlobalParams globalParams);

  void InitImageMosaicByBlocks(ModelParams input_img_params,
                               std::vector<ModelParams> overlap_img_params,
                               GlobalParams globalParams);

  void UpdateImageMosaic(ModelParams input_img_params,
                         std::vector<ModelParams> overlap_img_params,
                         GlobalParams globalParams);

  //albedo mosaic functions
  void InitAlbedoMosaic(ModelParams input_img_params,
                        std::vector<ModelParams> overlap_img_params,
                        GlobalParams globalParams);

  void UpdateAlbedoMosaic(ModelParams input_img_params,
                          std::vector<ModelParams> overlap_img_params,
                          GlobalParams globalParams);

  //reconstruction error functions
  void ComputeReconstructionErrorMap(ModelParams input_img_params,
                                     std::vector<ModelParams> overlap_img_params,
                                     GlobalParams globalParams,
                                     float *avgError, int *totalNumSamples);

}} // end vw::photometry

#endif//__VW_PHOTOMETRY_ALBEDO_H__
