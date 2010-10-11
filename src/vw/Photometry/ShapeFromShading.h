// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ShapeFromShading.h

#ifndef __VW_PHOTOMETRY_SHAPE_FROM_SHADING_H__
#define __VW_PHOTOMETRY_SHAPE_FROM_SHADING_H__

#include <vw/Image/ImageView.h>
#include <vw/Math/Vector.h>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {

  void UpdateHeightMap(ModelParams input_img_params, std::vector<ModelParams> overlap_img_params, GlobalParams globalParams);

  // Does conjugate gradient descent on the DEM, keeping all else fixed.
  void optimize_conjugate_gradient(ImageView<PixelGray<double> > *image_predicted,
                                   ImageView<PixelGray<double> > *image,
                                   ImageView<PixelGray<double> > *dem,
                                   ImageView<PixelGray<double> > *init_dem,
                                   ImageView<PixelGray<double> > *albedo,
                                   Vector3 *light_direction);

}} //end vw::photometry

#endif//__VW_PHOTOMETRY_SHAPE_FROM_SHADING_H__
